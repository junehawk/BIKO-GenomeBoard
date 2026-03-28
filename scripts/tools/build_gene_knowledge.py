#!/usr/bin/env python3
"""Build gene_knowledge.json from curated sources.

Source priority per gene:
1. CPIC (if PGx gene) → content_status: curated-cpic
2. CIViC local DB (if description exists) → curated-civic
3. NCBI Gene API (summary) → curated-ncbi
4. Minimal entry → auto-minimal

Additional enrichment (all genes):
- GeneReviews PMID → references
- ClinGen validity → finding_summary
- CIViC treatment evidence → treatment_strategies
"""

import json
import logging
from pathlib import Path
from typing import Dict, List, Optional
from scripts.common.api_utils import fetch_with_retry
from scripts.common.config import get
from scripts.db.query_civic import (
    get_gene_summary, get_treatment_summary, get_variant_evidence,
)
from scripts.tools.sources.ncbi_gene import fetch_gene_summary
from scripts.tools.sources.genreviews import fetch_genreviews_info
from scripts.db.query_orphanet import get_prevalence_text
from scripts.db.query_genreviews import get_genreviews_for_gene as get_genreviews_for_gene_local
from scripts.db.query_omim_mapping import get_mim_for_gene

logger = logging.getLogger(__name__)

CPIC_API = "https://api.cpicpgx.org/v1"
CPIC_PGX_GENES = [
    "CYP2D6", "CYP2C19", "CYP2C9", "CYP3A5", "DPYD",
    "NUDT15", "TPMT", "UGT1A1", "SLCO1B1", "VKORC1",
    "HLA-B", "HLA-A",
]


def fetch_cpic_gene(gene: str) -> Optional[Dict]:
    """Fetch gene info from CPIC API (unchanged from existing)."""
    data = fetch_with_retry(f"{CPIC_API}/gene?symbol=eq.{gene}&select=*")
    if not data:
        return None
    record = data[0] if isinstance(data, list) and data else data
    if not record or not record.get("symbol"):
        return None

    guidelines = fetch_with_retry(
        f"{CPIC_API}/pair?genesymbol=eq.{gene}&select=drugname,cpicLevel,pgkbLevel,guideline(name,url)"
    )
    drugs, refs = [], []
    if guidelines and isinstance(guidelines, list):
        for g in guidelines:
            drug = g.get("drugname", "")
            if drug:
                drugs.append(drug)
            gl = g.get("guideline")
            if gl and isinstance(gl, dict) and gl.get("name"):
                refs.append({"source": f"CPIC Guideline: {gl['name']}",
                             "note": f"CPIC Level {g.get('cpicLevel', 'N/A')}", "pmid": ""})

    treatment = f"CPIC guidelines available for: {', '.join(drugs)}." if drugs else ""
    return {
        "gene": gene, "full_name": record.get("name", ""),
        "function_summary": record.get("functionExampleSubstratesDrugs", ""),
        "clinical_significance": "Pharmacogenomically relevant gene with CPIC guidelines.",
        "associated_conditions": [f"{d} response" for d in drugs[:5]],
        "treatment_strategies": treatment, "frequency_prognosis": "",
        "finding_summary": f"{gene} is a pharmacogene with CPIC-level evidence for drug dosing.",
        "korean_specific_note": None, "hgvs": record.get("hgvs", {}),
        "references": refs or [{"source": "CPIC (cpicpgx.org)", "pmid": "", "note": "Auto-sourced"}],
        "content_status": "curated-cpic",
    }


def _try_clingen_validity(gene: str) -> Optional[str]:
    """Try local ClinGen DB for gene validity."""
    try:
        from scripts.db.query_local_clingen import get_gene_validity_local
        return get_gene_validity_local(gene)
    except Exception:
        return None


# Alias for monkeypatching in tests
get_gene_validity_local = _try_clingen_validity


def _build_gene_entry(gene: str) -> Dict:
    """Build a single gene knowledge entry using source priority chain."""

    # 1. CPIC (PGx genes)
    if gene in CPIC_PGX_GENES:
        cpic = fetch_cpic_gene(gene)
        if cpic:
            logger.info(f"  {gene}: CPIC (PGx)")
            return cpic

    # 2. CIViC local DB (cancer genes with descriptions)
    civic = get_gene_summary(gene)
    civic_description = civic.get("description", "") if civic else ""

    # 3. NCBI Gene API (universal fallback)
    ncbi = fetch_gene_summary(gene)

    # 4. GeneReviews PMID
    genreviews = fetch_genreviews_info(gene)

    # 5. ClinGen validity
    clingen = get_gene_validity_local(gene)

    # 6. CIViC treatment evidence
    treatment = get_treatment_summary(gene)

    # 7. CIViC evidence references
    civic_evidence = get_variant_evidence(gene)
    refs = []
    if civic_evidence:
        for e in civic_evidence[:5]:
            if e.get("pmid"):
                refs.append({"pmid": e["pmid"], "source": e.get("citation", ""),
                             "note": f"{e.get('evidence_type', '')} — {e.get('significance', '')}"})
    if genreviews:
        refs.append({"pmid": genreviews["pmid"], "source": "GeneReviews",
                      "note": genreviews.get("title", "")})

    # 8. Orphanet prevalence → frequency_prognosis
    orphanet_text = get_prevalence_text(gene)

    # 9. GeneReviews local DB (replaces API call when available)
    genreviews_local = get_genreviews_for_gene_local(gene)
    if genreviews_local and not genreviews:
        genreviews = genreviews_local
        refs.append({"pmid": genreviews_local["pmid"], "source": "GeneReviews",
                      "note": genreviews_local.get("title", "")})

    # 10. OMIM MIM mapping
    omim = get_mim_for_gene(gene)
    if omim:
        refs.append({"pmid": "", "source": f"OMIM #{omim['mim_number']}",
                      "note": omim.get("url", "")})

    # Compose finding_summary
    if civic_description:
        finding_summary = civic_description[:500]
        content_status = "curated-civic"
        source_label = "CIViC"
    elif ncbi and ncbi.get("summary"):
        finding_summary = ncbi["summary"][:500]
        content_status = "curated-ncbi"
        source_label = "NCBI Gene"
    else:
        finding_summary = f"{gene} — limited data available."
        content_status = "auto-minimal"
        source_label = "minimal"

    # Append ClinGen validity to finding_summary
    if clingen:
        finding_summary += f" ClinGen gene-disease validity: {clingen}."

    full_name = ""
    if ncbi:
        full_name = ncbi.get("full_name", "")
    elif civic:
        full_name = civic.get("aliases", "")

    logger.info(f"  {gene}: {source_label}" +
                (f" + GeneReviews" if genreviews else "") +
                (f" + ClinGen:{clingen}" if clingen else ""))

    return {
        "gene": gene,
        "full_name": full_name,
        "function_summary": ncbi.get("summary", "")[:300] if ncbi else "",
        "clinical_significance": "",
        "associated_conditions": [],
        "treatment_strategies": treatment or "",
        "frequency_prognosis": orphanet_text or "",
        "finding_summary": finding_summary,
        "korean_specific_note": None,
        "hgvs": {},
        "references": refs,
        "content_status": content_status,
    }


def build_knowledge(
    gene_list: List[str],
    output_path: str = "data/gene_knowledge.json",
) -> str:
    """Build gene_knowledge.json from curated sources."""
    logger.info(f"Building knowledge for {len(gene_list)} genes...")

    genes = []
    for gene in gene_list:
        entry = _build_gene_entry(gene)
        genes.append(entry)

    result = {"genes": genes}
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(result, f, indent=2, ensure_ascii=False)

    curated = sum(1 for g in genes if "curated" in g.get("content_status", ""))
    logger.info(f"Gene knowledge built: {len(genes)} genes ({curated} curated) → {output_path}")
    return output_path


def _load_gene_list_from_vcf(vcf_path: str) -> List[str]:
    """Extract unique gene symbols from a VCF file."""
    from scripts.intake.parse_vcf import parse_vcf
    variants = parse_vcf(vcf_path)
    return sorted(set(v.gene for v in variants if v.gene))


def _load_oncokb_genes() -> List[str]:
    """Load gene list from OncoKB cancer genes JSON."""
    path = get("paths.oncokb_genes", "data/oncokb_cancer_genes.json")
    try:
        with open(path) as f:
            data = json.load(f)
        return sorted(data.get("genes", {}).keys())
    except (FileNotFoundError, json.JSONDecodeError):
        return []


if __name__ == "__main__":
    import argparse
    logging.basicConfig(level=logging.INFO)
    parser = argparse.ArgumentParser(description="Build gene_knowledge.json from curated sources")
    parser.add_argument("--genes", help="Comma-separated gene list")
    parser.add_argument("--vcf", help="Extract genes from VCF file")
    parser.add_argument("--source", choices=["oncokb", "all"], help="Use predefined gene list")
    parser.add_argument("--output", default="data/gene_knowledge.json")
    args = parser.parse_args()

    if args.genes:
        gene_list = [g.strip() for g in args.genes.split(",")]
    elif args.vcf:
        gene_list = _load_gene_list_from_vcf(args.vcf)
    elif args.source == "oncokb":
        gene_list = _load_oncokb_genes()
    elif args.source == "all":
        gene_list = _load_oncokb_genes() + CPIC_PGX_GENES
        gene_list = sorted(set(gene_list))
    else:
        parser.error("Specify --genes, --vcf, or --source")

    build_knowledge(gene_list, args.output)
