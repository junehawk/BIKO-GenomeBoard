"""OMIM gene-disease association lookup."""
import logging
from typing import Dict, List, Optional
from scripts.common.config import get

logger = logging.getLogger(__name__)

# Static OMIM data for common rare disease genes (avoid API key requirement)
# This will be expanded when OMIM API key is configured
OMIM_DATA = {
    "TP53": {"mim": "191170", "phenotypes": ["Li-Fraumeni syndrome"], "inheritance": "AD"},
    "BRCA2": {"mim": "600185", "phenotypes": ["Hereditary breast-ovarian cancer syndrome", "Fanconi anemia D1"], "inheritance": "AD/AR"},
    "CFTR": {"mim": "602421", "phenotypes": ["Cystic fibrosis", "CBAVD"], "inheritance": "AR"},
    "ATM": {"mim": "607585", "phenotypes": ["Ataxia-telangiectasia", "Breast cancer susceptibility"], "inheritance": "AR/AD"},
    "MUTYH": {"mim": "604933", "phenotypes": ["MUTYH-associated polyposis", "Colorectal adenomas"], "inheritance": "AR"},
    "PALB2": {"mim": "610355", "phenotypes": ["Hereditary breast cancer", "Fanconi anemia N"], "inheritance": "AD/AR"},
    "PTPN11": {"mim": "176876", "phenotypes": ["Noonan syndrome 1", "LEOPARD syndrome", "Juvenile myelomonocytic leukemia"], "inheritance": "AD"},
    "APOE": {"mim": "107741", "phenotypes": ["Alzheimer disease 2", "Hyperlipoproteinemia type III"], "inheritance": "AD/complex"},
    "NUDT15": {"mim": "615792", "phenotypes": ["Thiopurine sensitivity"], "inheritance": "AR"},
    "CYP2C19": {"mim": "124020", "phenotypes": ["Poor drug metabolism (CYP2C19)"], "inheritance": "AR"},
    "HLA-B": {"mim": "142830", "phenotypes": ["Drug hypersensitivity"], "inheritance": "codominant"},
}


def query_omim(gene: str) -> Optional[Dict]:
    """Get OMIM data for a gene. Returns dict with mim, phenotypes, inheritance."""
    return OMIM_DATA.get(gene)
