"""Clinical Board runner — orchestrates the full diagnostic synthesis process."""

import json
import logging
import os
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Optional, Union

from scripts.clinical_board.models import (
    AgentOpinion,
    BoardOpinion,
    CancerBoardOpinion,
)
from scripts.clinical_board.ollama_client import OllamaClient
from scripts.clinical_board.case_briefing import build_case_briefing
from scripts.clinical_board.domain_sheets import build_domain_sheet
from scripts.clinical_board.kb_query import query_prior_knowledge
from scripts.clinical_board.knowledge_base import KnowledgeBase
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
        from scripts.clinical_board.agents.therapeutic_target import TherapeuticTargetAnalyst
        from scripts.clinical_board.agents.tumor_genomics import TumorGenomicsSpecialist
        from scripts.clinical_board.agents.pgx_specialist import PGxSpecialist
        from scripts.clinical_board.agents.clinical_evidence import ClinicalEvidenceAnalyst
        return [
            TherapeuticTargetAnalyst(client=client, model=model, language=language),
            TumorGenomicsSpecialist(client=client, model=model, language=language),
            PGxSpecialist(client=client, model=model, language=language),
            ClinicalEvidenceAnalyst(client=client, model=model, language=language),
        ]

    from scripts.clinical_board.agents.variant_pathologist import VariantPathologist
    from scripts.clinical_board.agents.disease_geneticist import DiseaseGeneticist
    from scripts.clinical_board.agents.pgx_specialist import PGxSpecialist
    from scripts.clinical_board.agents.literature_analyst import LiteratureAnalyst
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
            logger.warning(f"[Clinical Board] Model '{model_name}' not found. "
                          f"Run 'ollama pull {model_name}' to download.")

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
            domain_sheet = build_domain_sheet(
                agent.domain, mode, report_data.get("variants", []), report_data
            )
            opinion = agent.analyze(
                briefing, domain_sheet=domain_sheet, prior_knowledge=prior_knowledge
            )
            opinions.append(opinion)
            logger.info(f"  [Board] {agent.agent_name}: {opinion.confidence} confidence, "
                       f"{len(opinion.findings)} findings")
        except Exception as e:
            logger.error(f"  [Board] {agent.agent_name} failed: {e}")
            opinions.append(AgentOpinion(
                agent_name=agent.agent_name, domain="error",
                findings=[{"finding": f"Agent failed: {e}", "evidence": "", "confidence": "low"}],
                confidence="low",
            ))

    # Step 3: Board Chair synthesis
    logger.info("[Clinical Board] Board Chair synthesizing opinions...")
    chair = _load_chair(client, chair_model, language)
    board_opinion = chair.synthesize(briefing, opinions, mode=mode)
    board_opinion.agent_opinions = opinions

    # Step 4: Save decisions to Knowledge Base
    if get("knowledge_base.enabled", False):
        try:
            kb = KnowledgeBase(kb_db, kb_path)
            saved_genes: set[str] = set()
            skipped = 0
            for v in report_data.get("variants", [])[:20]:
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
                    logger.warning(
                        f"[Clinical Board] KB save skipped {gene}/{variant_id}: {row_err}"
                    )
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
        f"[Clinical Board] Complete in {elapsed:.1f}s — "
        f"{label}: {summary} ({board_opinion.confidence} confidence)"
    )

    return board_opinion
