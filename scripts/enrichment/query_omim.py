"""OMIM gene-disease association lookup.

Uses genemap2.txt SQLite DB when available, falls back to static dict.
"""

import logging
from typing import Dict, Optional

logger = logging.getLogger(__name__)

# Static OMIM data fallback for common genes (used when genemap DB not available)
_STATIC_OMIM = {
    "TP53": {"mim": "191170", "phenotypes": ["Li-Fraumeni syndrome"], "inheritance": "AD"},
    "BRCA2": {
        "mim": "600185",
        "phenotypes": ["Hereditary breast-ovarian cancer syndrome", "Fanconi anemia D1"],
        "inheritance": "AD/AR",
    },
    "CFTR": {"mim": "602421", "phenotypes": ["Cystic fibrosis", "CBAVD"], "inheritance": "AR"},
    "ATM": {
        "mim": "607585",
        "phenotypes": ["Ataxia-telangiectasia", "Breast cancer susceptibility"],
        "inheritance": "AR/AD",
    },
    "MUTYH": {
        "mim": "604933",
        "phenotypes": ["MUTYH-associated polyposis", "Colorectal adenomas"],
        "inheritance": "AR",
    },
    "PALB2": {"mim": "610355", "phenotypes": ["Hereditary breast cancer", "Fanconi anemia N"], "inheritance": "AD/AR"},
    "PTPN11": {
        "mim": "176876",
        "phenotypes": ["Noonan syndrome 1", "LEOPARD syndrome", "Juvenile myelomonocytic leukemia"],
        "inheritance": "AD",
    },
    "APOE": {
        "mim": "107741",
        "phenotypes": ["Alzheimer disease 2", "Hyperlipoproteinemia type III"],
        "inheritance": "AD/complex",
    },
    "NUDT15": {"mim": "615792", "phenotypes": ["Thiopurine sensitivity"], "inheritance": "AR"},
    "CYP2C19": {"mim": "124020", "phenotypes": ["Poor drug metabolism (CYP2C19)"], "inheritance": "AR"},
    "HLA-B": {"mim": "142830", "phenotypes": ["Drug hypersensitivity"], "inheritance": "codominant"},
}


def query_omim(gene: str) -> Optional[Dict]:
    """Get OMIM data for a gene. Uses genemap DB first, then static fallback."""
    # Try genemap2 DB first
    try:
        from scripts.storage.query_omim_genemap import get_gene_phenotypes, get_inheritance_patterns

        phenotypes_data = get_gene_phenotypes(gene)
        if phenotypes_data:
            inheritance_list = get_inheritance_patterns(gene) or []
            return {
                "mim": phenotypes_data[0].get("mim_number", ""),
                "phenotypes": [p["phenotype"] for p in phenotypes_data if p.get("phenotype")],
                "inheritance": "/".join(inheritance_list) if inheritance_list else "",
                "source": "genemap2",
            }
    except ImportError:
        pass
    except Exception as e:
        logger.debug(f"genemap DB lookup failed for {gene}: {e}")

    # Fallback to static dict
    result = _STATIC_OMIM.get(gene)
    if result:
        return {**result, "source": "static"}
    return None
