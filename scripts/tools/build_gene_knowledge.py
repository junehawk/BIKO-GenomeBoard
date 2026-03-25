#!/usr/bin/env python3
"""Build gene_knowledge.json from curated sources (CPIC API + PubMed)."""

import json
import logging
from pathlib import Path
from typing import Dict, List, Optional
from scripts.common.api_utils import fetch_with_retry

logger = logging.getLogger(__name__)

CPIC_API = "https://api.cpicpgx.org/v1"
CPIC_PGX_GENES = [
    "CYP2D6", "CYP2C19", "CYP2C9", "CYP3A5", "DPYD",
    "NUDT15", "TPMT", "UGT1A1", "SLCO1B1", "VKORC1",
    "HLA-B", "HLA-A",
]


def fetch_cpic_gene(gene: str) -> Optional[Dict]:
    """Fetch gene info from CPIC API and format for gene_knowledge.json."""
    data = fetch_with_retry(f"{CPIC_API}/gene?symbol=eq.{gene}&select=*")
    if not data:
        return None

    # CPIC returns array
    record = data[0] if isinstance(data, list) and data else data
    if not record or not record.get("symbol"):
        return None

    # Fetch associated drug-gene pairs for guidelines
    guidelines = fetch_with_retry(
        f"{CPIC_API}/pair?genesymbol=eq.{gene}&select=drugname,cpicLevel,pgkbLevel,guideline(name,url)"
    )
    drugs = []
    refs = []
    if guidelines and isinstance(guidelines, list):
        for g in guidelines:
            drug = g.get("drugname", "")
            if drug:
                drugs.append(drug)
            gl = g.get("guideline")
            if gl and isinstance(gl, dict) and gl.get("name"):
                refs.append({
                    "source": f"CPIC Guideline: {gl['name']}",
                    "note": f"CPIC Level {g.get('cpicLevel', 'N/A')}",
                    "pmid": "",
                })

    treatment = f"CPIC guidelines available for: {', '.join(drugs)}." if drugs else ""

    return {
        "gene": gene,
        "full_name": record.get("name", ""),
        "function_summary": record.get("functionExampleSubstratesDrugs", ""),
        "clinical_significance": "Pharmacogenomically relevant gene with CPIC guidelines.",
        "associated_conditions": [f"{d} response" for d in drugs[:5]],
        "treatment_strategies": treatment,
        "frequency_prognosis": "",
        "finding_summary": f"{gene} is a pharmacogene with CPIC-level evidence for drug dosing.",
        "korean_specific_note": None,
        "hgvs": record.get("hgvs", {}),
        "references": refs or [{"source": "CPIC (cpicpgx.org)", "pmid": "", "note": "Auto-sourced"}],
        "content_status": "curated-cpic",
    }


def build_knowledge(
    existing_path: str,
    output_path: str,
    cpic_genes: Optional[List[str]] = None,
) -> str:
    """Merge CPIC-curated data into existing gene_knowledge.json."""
    cpic_genes = cpic_genes or CPIC_PGX_GENES

    # Load existing
    existing = {"genes": []}
    if Path(existing_path).exists():
        with open(existing_path) as f:
            existing = json.load(f)

    gene_map = {g["gene"]: g for g in existing["genes"]}

    # Fetch CPIC data for PGx genes
    for gene in cpic_genes:
        logger.info(f"Fetching CPIC data for {gene}...")
        cpic_data = fetch_cpic_gene(gene)
        if cpic_data:
            if gene in gene_map:
                # Preserve fields CPIC doesn't provide
                for key in ["korean_specific_note", "hgvs"]:
                    if not cpic_data.get(key) and gene_map[gene].get(key):
                        cpic_data[key] = gene_map[gene][key]
            gene_map[gene] = cpic_data
            logger.info(f"  → {gene}: updated from CPIC")
        else:
            logger.warning(f"  → {gene}: CPIC data not available, keeping existing")

    result = {"genes": list(gene_map.values())}
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(result, f, indent=2, ensure_ascii=False)

    curated = sum(1 for g in result["genes"] if "curated" in g.get("content_status", ""))
    logger.info(f"Gene knowledge built: {len(result['genes'])} genes ({curated} curated) → {output_path}")
    return output_path


if __name__ == "__main__":
    import argparse
    logging.basicConfig(level=logging.INFO)
    parser = argparse.ArgumentParser()
    parser.add_argument("--existing", default="data/gene_knowledge.json")
    parser.add_argument("--output", default="data/gene_knowledge.json")
    args = parser.parse_args()
    build_knowledge(args.existing, args.output)
