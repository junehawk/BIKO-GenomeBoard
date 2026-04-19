"""Clinical Board runner — orchestrates the full diagnostic synthesis process."""

import json
import logging
import os
import time
from typing import Optional, Union

from scripts.clinical_board.case_briefing import build_case_briefing
from scripts.clinical_board.curated_treatments import curate_treatments
from scripts.clinical_board.domain_sheets import build_domain_sheet
from scripts.clinical_board.kb_query import query_prior_knowledge
from scripts.clinical_board.knowledge_base import KnowledgeBase
from scripts.clinical_board.models import (
    AgentOpinion,
    BoardOpinion,
    CancerBoardOpinion,
)
from scripts.clinical_board.narrative_scrubber import scrub_opinion
from scripts.clinical_board.ollama_client import OllamaClient
from scripts.clinical_board.template_renderer_chair import render_from_curated
from scripts.clinical_board.variant_selector import select_board_variants
from scripts.common.config import get

logger = logging.getLogger(__name__)


def _load_agents(
    client: OllamaClient,
    model: str,
    language: str,
    mode: str = "rare-disease",
):
    """Lazy-load domain agents for the given mode."""
    if mode == "cancer":
        from scripts.clinical_board.agents.clinical_evidence import ClinicalEvidenceAnalyst
        from scripts.clinical_board.agents.pgx_specialist import PGxSpecialist
        from scripts.clinical_board.agents.therapeutic_target import TherapeuticTargetAnalyst
        from scripts.clinical_board.agents.tumor_genomics import TumorGenomicsSpecialist

        return [
            TherapeuticTargetAnalyst(client=client, model=model, language=language),
            TumorGenomicsSpecialist(client=client, model=model, language=language),
            PGxSpecialist(client=client, model=model, language=language),
            ClinicalEvidenceAnalyst(client=client, model=model, language=language),
        ]

    from scripts.clinical_board.agents.disease_geneticist import DiseaseGeneticist
    from scripts.clinical_board.agents.literature_analyst import LiteratureAnalyst
    from scripts.clinical_board.agents.pgx_specialist import PGxSpecialist
    from scripts.clinical_board.agents.variant_pathologist import VariantPathologist

    return [
        VariantPathologist(client=client, model=model, language=language),
        DiseaseGeneticist(client=client, model=model, language=language),
        PGxSpecialist(client=client, model=model, language=language),
        LiteratureAnalyst(client=client, model=model, language=language),
    ]


def _load_chair(client: OllamaClient, model: str, language: str):
    from scripts.clinical_board.agents.board_chair import BoardChair

    return BoardChair(client=client, model=model, language=language)


def _summarize_context(report_data: dict) -> str:
    """Compact summary of clinical context for KB storage."""
    parts = []
    sample_id = report_data.get("sample_id")
    if sample_id:
        parts.append(f"sample:{sample_id}")
    note = report_data.get("clinical_note", "")
    if note:
        snippet = note.strip().replace("\n", " ")[:120]
        parts.append(f"note:{snippet}")
    total = report_data.get("summary", {}).get("total")
    if total is not None:
        parts.append(f"variants:{total}")
    return "|".join(parts)


def run_clinical_board(
    report_data: dict,
    mode: str,
    ollama_url: str = None,
    agent_model: str = None,
    chair_model: str = None,
    language: str = None,
) -> Optional[Union[BoardOpinion, CancerBoardOpinion]]:
    """Run the full Clinical Board analysis.

    Args:
        report_data: Output from run_pipeline()
        mode: "cancer" or "rare-disease"
        ollama_url: Ollama server URL (default from config)
        agent_model: Model for domain agents (default: medgemma:27b)
        chair_model: Model for Board Chair (default: gemma4:31b)

    Returns:
        BoardOpinion with synthesized diagnostic opinion, or None on failure.
    """
    start = time.time()
    agent_model = agent_model or get("clinical_board.agent_model", "medgemma:27b")
    chair_model = chair_model or get("clinical_board.chair_model", "alibayram/medgemma:27b")
    language = language or get("clinical_board.language", "en")

    # Initialize client
    client = OllamaClient(base_url=ollama_url)
    if not client.is_available():
        logger.error("[Clinical Board] Ollama server is not available. Skipping.")
        return None

    # Check models
    for model_name in [agent_model, chair_model]:
        if not client.has_model(model_name):
            logger.warning(
                f"[Clinical Board] Model '{model_name}' not found. Run 'ollama pull {model_name}' to download."
            )

    # Step 0: Pre-filter variants with the deterministic selector. Every
    # downstream consumer (case briefing, domain sheets, KB save, render
    # footer) reads the same filtered list via report_data["_board_variants"]
    # so the audit trail stays consistent.
    raw_variants = report_data.get("variants", []) or []
    board_variants, selection_metadata = select_board_variants(raw_variants, mode, report_data=report_data)
    report_data["_board_variants"] = board_variants
    report_data["_board_selection_metadata"] = selection_metadata

    # Write selection_reason_list back to report_data["variants"] so the
    # rare-disease template's denovo_badge macro can render it. The
    # selector runs AFTER build_variant_records so mutating the Variant
    # object is not enough — we need to update the already-snapshotted
    # variant record dicts. Match by variant_id.
    _reason_map = {
        (bv.get("variant") or bv.get("variant_id") or ""): bv.get("selection_reason_list", [])
        for bv in board_variants
        if bv.get("selection_reason_list")
    }
    if _reason_map:
        for v_record in raw_variants:
            vid = v_record.get("variant") or v_record.get("variant_id") or ""
            if vid in _reason_map:
                v_record["selection_reason_list"] = _reason_map[vid]

    # Promote board-admitted VUS variants to detailed_variants. The
    # classify.py split_variants_for_display() runs BEFORE the Board, so
    # a VUS admitted by the selector (and potentially cited by the Chair)
    # ends up in the omitted_variants appendix instead of the main detail
    # pages. Observed on ASD-10293: CHD8 VUS was the Chair's primary
    # diagnosis ("CHD8-related neurodevelopmental disorder") but the
    # per-variant detail pages did not include CHD8 at all — the reviewer
    # saw the diagnosis with no variant context to back it up. Moving
    # board-admitted variants into detailed restores that context. The
    # move is only from omitted → detailed (never the reverse), so P/LP
    # already in detailed stay untouched.
    detailed = report_data.get("detailed_variants") or []
    omitted = report_data.get("omitted_variants") or []
    board_admitted_ids = {
        (bv.get("variant") or bv.get("variant_id") or "")
        for bv in board_variants
        if (bv.get("variant") or bv.get("variant_id"))
    }
    if board_admitted_ids and detailed is not None and omitted:
        detailed_ids = {(v.get("variant") or v.get("variant_id") or "") for v in detailed}
        to_promote = [
            v
            for v in omitted
            if (v.get("variant") or v.get("variant_id") or "") in board_admitted_ids
            and (v.get("variant") or v.get("variant_id") or "") not in detailed_ids
        ]
        if to_promote:
            report_data["detailed_variants"] = list(detailed) + to_promote
            report_data["omitted_variants"] = [v for v in omitted if v not in to_promote]
            logger.info(
                "[Clinical Board] Promoted %d board-admitted VUS variant(s) from omitted to detailed pages",
                len(to_promote),
            )

    # Re-sort report_data["variants"] so the summary-page short-lists
    # (Candidate Gene Ranking top-15, VUS panel top-5, No Reportable
    # top-10) surface clinically interesting variants first. Previously
    # the order was genomic coordinate, so a Chair-cited variant on
    # chr14 (CHD8 was position #36) would fall below 35 random chr1/2
    # intergenic VUS and never appear in any page-1 section.
    #
    # Sort key (ascending; lower = higher priority):
    #   1. Board-admitted first (selection_reason_list non-empty)
    #   2. Classification rank (P < LP < Drug Response < VUS < LB < B)
    #   3. HPO score (negated — higher HPO first)
    #   4. Has gene annotation (gene-bearing before gene=None)
    #   5. variant_id (stable genomic-coordinate tiebreaker)
    _CLASS_RANK = {
        "pathogenic": 0,
        "likely pathogenic": 1,
        "drug response": 2,
        "risk factor": 3,
        "vus": 4,
        "likely benign": 5,
        "benign": 6,
    }

    def _variant_sort_key(v: dict) -> tuple:
        board_admitted = 0 if v.get("selection_reason_list") else 1
        cls_rank = _CLASS_RANK.get((v.get("classification") or "VUS").lower(), 9)
        hpo_score = -(v.get("hpo_score") or 0)
        has_gene = 0 if v.get("gene") else 1
        return (board_admitted, cls_rank, hpo_score, has_gene, v.get("variant") or "")

    if raw_variants:
        raw_variants.sort(key=_variant_sort_key)

    logger.info(
        f"[Clinical Board] Pre-filter: {selection_metadata['total_input']} → "
        f"{selection_metadata['selected']} variants "
        f"({selection_metadata['criteria_summary']})"
    )

    # v2.2 · Curate-then-narrate — deterministic treatment resolver.
    # Cancer mode only; rare-disease bypasses the OncoKB/CIViC curator.
    if mode == "cancer":
        try:
            curated = curate_treatments(
                board_variants,
                offline_mode=get("clinical_board.curated_treatments.offline_mode", False),
            )
            report_data["_curated_treatments"] = curated
            total_rows = sum(len(rows) for rows in curated.values())
            logger.info(f"[Clinical Board] Curated treatments: {len(curated)} variants → {total_rows} rows")
        except Exception as curator_err:
            logger.warning(f"[Clinical Board] Curator failed: {curator_err}")
            report_data["_curated_treatments"] = {}

    # Step 1: Build case briefing
    logger.info("[Clinical Board] Building case briefing...")
    briefing = build_case_briefing(report_data, mode)
    logger.debug(f"[Clinical Board] Briefing length: {len(briefing)} chars")

    # Query KB for prior knowledge
    variant_ids = [v.get("variant", "") for v in report_data.get("variants", [])]
    kb_path = get("knowledge_base.path", "data/knowledge_base")
    kb_db = os.path.join(kb_path, "kb.sqlite3")
    prior_knowledge = ""
    if get("knowledge_base.enabled", False) and os.path.exists(kb_db):
        prior_knowledge = query_prior_knowledge(kb_db, variant_ids, mode)

    # Step 2: Run domain agents in parallel
    logger.info("[Clinical Board] Running domain agents...")
    agents = _load_agents(client, agent_model, language, mode=mode)
    opinions: list[AgentOpinion] = []

    # Run agents sequentially on single GPU to avoid resource contention
    # (parallel execution causes timeouts when agents compete for GPU)
    for agent in agents:
        try:
            domain_sheet = build_domain_sheet(agent.domain, mode, board_variants, report_data)
            opinion = agent.analyze(briefing, domain_sheet=domain_sheet, prior_knowledge=prior_knowledge)
            opinions.append(opinion)
            logger.info(
                f"  [Board] {agent.agent_name}: {opinion.confidence} confidence, {len(opinion.findings)} findings"
            )
        except Exception as e:
            logger.error(f"  [Board] {agent.agent_name} failed: {e}")
            opinions.append(
                AgentOpinion(
                    agent_name=agent.agent_name,
                    domain="error",
                    findings=[{"finding": f"Agent failed: {e}", "evidence": "", "confidence": "low"}],
                    confidence="low",
                )
            )

    # Step 3: Board Chair synthesis
    logger.info("[Clinical Board] Board Chair synthesizing opinions...")
    chair = _load_chair(client, chair_model, language)
    curated_for_chair = report_data.get("_curated_treatments", {}) if mode == "cancer" else None
    board_opinion = chair.synthesize(briefing, opinions, curated_treatments=curated_for_chair, mode=mode)
    board_opinion.agent_opinions = opinions
    # Propagate selection audit trail for render.py footer and KB traceability.
    board_opinion.selection_metadata = selection_metadata

    # Step 3b (cancer only): narrative_scrubber patient-safety gate.
    # Drops any (curated_id, variant_key) pair the curator did not emit and
    # strips banned drug tokens out of every prose field.
    #
    # The template_renderer_chair fallback fires when, after scrubbing, the
    # opinion has zero treatment rows AND the curator did emit non-empty
    # rows. Two separate failure modes both hit this condition:
    #   (a) the LLM produced zero rows — ignored the curated block entirely.
    #   (b) the LLM produced rows but every row was dropped by the scrubber
    #       (invalid curated_id / variant_key pair — schema drift on small
    #       local models, or cross-variant paste attack).
    # In both cases the fallback renderer is safer than showing an empty
    # therapy table because it surfaces the curator's authoritative rows
    # with an explicit "fallback / research reference library" label and
    # confidence=low. The fallback never emits a row the LLM endorsed; it
    # only surfaces what the deterministic curator already produced.
    if mode == "cancer":
        pre_scrub_count = len(getattr(board_opinion, "treatment_options", []) or [])
        curated_nonempty = bool(curated_for_chair and any(rows for rows in curated_for_chair.values()))
        try:
            stats = scrub_opinion(board_opinion, curated_for_chair or {})
            logger.info(
                f"[Clinical Board] narrative_scrubber: kept={stats['kept']} "
                f"dropped={stats['dropped']} banned_terms={stats['banned_terms']}"
            )
        except Exception as scrub_err:
            logger.warning(f"[Clinical Board] narrative_scrubber failed: {scrub_err}")

        post_scrub_count = len(getattr(board_opinion, "treatment_options", []) or [])
        if curated_nonempty and post_scrub_count == 0:
            if pre_scrub_count == 0:
                reason = "LLM produced zero treatment rows"
            else:
                reason = (
                    f"LLM produced {pre_scrub_count} rows but scrubber dropped "
                    f"all of them (schema drift or paste attack)"
                )
            logger.warning(
                "[Clinical Board] Board Chair fallback triggered — %s. "
                "Filling treatment rows from deterministic curator via "
                "template_renderer_chair (hybrid merge — Chair narrative preserved).",
                reason,
            )
            fallback = render_from_curated(curated_for_chair or {}, agent_opinions=opinions)

            # Hybrid merge: the Chair's narrative synthesis (therapeutic
            # implications, actionable findings, clinical actions, monitoring
            # plan, immunotherapy assessment) is the LLM's primary value-add
            # over the raw agent opinions. Throwing it all away because
            # treatment_options was empty/invalid wastes ~4 min of GPU time
            # and produces a generic "unable to produce recommendations"
            # placeholder instead of actual clinical reasoning.
            #
            # We keep the Chair's narrative and fill ONLY the treatment rows
            # from the deterministic fallback. The scrubber already ran on
            # the narrative fields so hallucinated drug names (if any) are
            # redacted. Set consensus to "hybrid-fallback" so downstream
            # render.py and the reviewer can tell which parts are LLM
            # (narrative) vs deterministic (treatment rows).
            chair_has_narrative = bool(
                getattr(board_opinion, "therapeutic_implications", "")
                and len(getattr(board_opinion, "therapeutic_implications", "")) > 30
            )
            if chair_has_narrative:
                board_opinion.treatment_options = fallback.treatment_options
                board_opinion.agent_consensus = "hybrid-fallback"
                board_opinion.selection_metadata = selection_metadata
            else:
                # Chair produced no useful narrative either — full fallback.
                fallback.selection_metadata = selection_metadata
                board_opinion = fallback

    # Step 4: Save decisions to Knowledge Base
    if get("knowledge_base.enabled", False):
        try:
            kb = KnowledgeBase(kb_db, kb_path)
            saved_genes: set[str] = set()
            skipped = 0
            for v in board_variants:
                gene = v.get("gene") or ""
                variant_id = v.get("variant") or ""
                if not gene or not variant_id:
                    skipped += 1
                    continue
                try:
                    kb.save_decision(
                        sample_id=report_data.get("sample_id", "unknown"),
                        mode=mode,
                        gene=gene,
                        variant=variant_id,
                        hgvsp=v.get("hgvsp", ""),
                        classification=v.get("classification", ""),
                        board_diagnosis=getattr(
                            board_opinion,
                            "primary_diagnosis",
                            getattr(board_opinion, "therapeutic_implications", ""),
                        ),
                        board_confidence=board_opinion.confidence,
                        clinical_context_summary=_summarize_context(report_data),
                        agent_consensus=board_opinion.agent_consensus,
                        raw_opinion_json=json.dumps(board_opinion.__dict__, default=str),
                    )
                    saved_genes.add(gene)
                except Exception as row_err:
                    logger.warning(f"[Clinical Board] KB save skipped {gene}/{variant_id}: {row_err}")
                    skipped += 1
            for gene in saved_genes:
                kb.generate_gene_wiki(gene, mode)
            if saved_genes:
                kb.update_log(report_data.get("sample_id"), mode)
            if skipped:
                logger.info(f"[Clinical Board] KB: {len(saved_genes)} genes saved, {skipped} rows skipped")
        except Exception as e:
            logger.warning(f"[Clinical Board] KB save failed: {e}")

    elapsed = time.time() - start
    if isinstance(board_opinion, CancerBoardOpinion):
        summary = board_opinion.therapeutic_implications
        label = "Therapeutic implications"
    else:
        summary = board_opinion.primary_diagnosis
        label = "Primary diagnosis"
    logger.info(
        f"[Clinical Board] Complete in {elapsed:.1f}s — {label}: {summary} ({board_opinion.confidence} confidence)"
    )

    return board_opinion
