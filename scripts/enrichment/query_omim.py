"""OMIM gene-disease association lookup.

Uses genemap2.txt SQLite DB when available, falls back to static dict.
"""

import logging
from typing import Dict, Optional

logger = logging.getLogger(__name__)

# Static OMIM data fallback (used when genemap2.txt SQLite DB is absent).
#
# v2.6 (2026-04-30): expanded from 11 to ~80 genes. Coverage is anchored on
# the ACMG SF v3.2 secondary-findings gene list plus high-yield Lynch /
# HBOC / polyposis / neurodevelopmental / cardiomyopathy / connective-tissue
# / RASopathy / hematology / pharmacogenomics genes seen in real Korean
# WGS workloads. Each entry was cross-checked against OMIM gene MIM
# numbers (public catalogue, no licensed text reused). Phenotype labels
# are short, ASCII-only summaries; richer free-text annotation requires
# the licensed genemap2.txt build.
#
# When OMIM access lands, the SQLite-backed lookup
# (``scripts.storage.query_omim_genemap``) supersedes this dict. Until
# then, the static fallback gives the rare-disease report a meaningful
# ``inheritance`` column on most actionable genes instead of leaving it
# blank.
_STATIC_OMIM = {
    # --- Hereditary breast / ovarian cancer (HBOC) -------------------------
    "BRCA1": {"mim": "113705", "phenotypes": ["Hereditary breast-ovarian cancer 1"], "inheritance": "AD"},
    "BRCA2": {
        "mim": "600185",
        "phenotypes": ["Hereditary breast-ovarian cancer 2", "Fanconi anemia D1"],
        "inheritance": "AD/AR",
    },
    "PALB2": {"mim": "610355", "phenotypes": ["Hereditary breast cancer", "Fanconi anemia N"], "inheritance": "AD/AR"},
    "ATM": {
        "mim": "607585",
        "phenotypes": ["Ataxia-telangiectasia", "Breast cancer susceptibility"],
        "inheritance": "AR/AD",
    },
    "CHEK2": {"mim": "604373", "phenotypes": ["Breast cancer susceptibility"], "inheritance": "AD"},
    "BARD1": {"mim": "601593", "phenotypes": ["Breast cancer susceptibility"], "inheritance": "AD"},
    "BRIP1": {
        "mim": "605882",
        "phenotypes": ["Fanconi anemia J", "Breast cancer susceptibility"],
        "inheritance": "AD/AR",
    },
    "RAD51C": {
        "mim": "602774",
        "phenotypes": ["Fanconi anemia O", "Breast-ovarian cancer susceptibility"],
        "inheritance": "AD/AR",
    },
    "RAD51D": {"mim": "602954", "phenotypes": ["Breast-ovarian cancer susceptibility"], "inheritance": "AD"},
    # --- Lynch syndrome / mismatch repair ---------------------------------
    "MLH1": {"mim": "120436", "phenotypes": ["Lynch syndrome"], "inheritance": "AD"},
    "MSH2": {"mim": "609309", "phenotypes": ["Lynch syndrome"], "inheritance": "AD"},
    "MSH6": {"mim": "600678", "phenotypes": ["Lynch syndrome"], "inheritance": "AD"},
    "PMS2": {
        "mim": "600259",
        "phenotypes": ["Lynch syndrome", "Constitutional mismatch repair deficiency"],
        "inheritance": "AD/AR",
    },
    "EPCAM": {"mim": "185535", "phenotypes": ["Lynch syndrome (EPCAM deletions)"], "inheritance": "AD"},
    # --- Colorectal polyposis ---------------------------------------------
    "APC": {"mim": "611731", "phenotypes": ["Familial adenomatous polyposis"], "inheritance": "AD"},
    "MUTYH": {
        "mim": "604933",
        "phenotypes": ["MUTYH-associated polyposis", "Colorectal adenomas"],
        "inheritance": "AR",
    },
    "POLE": {"mim": "174762", "phenotypes": ["Polymerase proofreading-associated polyposis"], "inheritance": "AD"},
    "POLD1": {"mim": "174761", "phenotypes": ["Polymerase proofreading-associated polyposis"], "inheritance": "AD"},
    "STK11": {"mim": "602216", "phenotypes": ["Peutz-Jeghers syndrome"], "inheritance": "AD"},
    "SMAD4": {
        "mim": "600993",
        "phenotypes": ["Juvenile polyposis", "Hereditary hemorrhagic telangiectasia"],
        "inheritance": "AD",
    },
    "BMPR1A": {"mim": "601299", "phenotypes": ["Juvenile polyposis"], "inheritance": "AD"},
    # --- Other tumor-predisposition syndromes -----------------------------
    "TP53": {"mim": "191170", "phenotypes": ["Li-Fraumeni syndrome"], "inheritance": "AD"},
    "CDH1": {
        "mim": "192090",
        "phenotypes": ["Hereditary diffuse gastric cancer", "Lobular breast cancer"],
        "inheritance": "AD",
    },
    "VHL": {"mim": "608537", "phenotypes": ["von Hippel-Lindau syndrome"], "inheritance": "AD"},
    "RET": {
        "mim": "164761",
        "phenotypes": ["Multiple endocrine neoplasia 2A/2B", "Hirschsprung disease"],
        "inheritance": "AD",
    },
    "MEN1": {"mim": "613733", "phenotypes": ["Multiple endocrine neoplasia 1"], "inheritance": "AD"},
    "NF1": {"mim": "613113", "phenotypes": ["Neurofibromatosis 1"], "inheritance": "AD"},
    "NF2": {"mim": "607379", "phenotypes": ["Neurofibromatosis 2"], "inheritance": "AD"},
    "PTEN": {"mim": "601728", "phenotypes": ["Cowden syndrome", "PTEN hamartoma tumor syndrome"], "inheritance": "AD"},
    "RB1": {"mim": "614041", "phenotypes": ["Retinoblastoma"], "inheritance": "AD"},
    "TSC1": {"mim": "605284", "phenotypes": ["Tuberous sclerosis 1"], "inheritance": "AD"},
    "TSC2": {"mim": "191092", "phenotypes": ["Tuberous sclerosis 2"], "inheritance": "AD"},
    # --- Cardiomyopathy / arrhythmia (ACMG SF v3.2 core) ------------------
    "MYH7": {
        "mim": "160760",
        "phenotypes": ["Hypertrophic cardiomyopathy", "Dilated cardiomyopathy"],
        "inheritance": "AD",
    },
    "MYBPC3": {"mim": "600958", "phenotypes": ["Hypertrophic cardiomyopathy"], "inheritance": "AD"},
    "TNNT2": {"mim": "191045", "phenotypes": ["Hypertrophic / dilated cardiomyopathy"], "inheritance": "AD"},
    "TNNI3": {"mim": "191044", "phenotypes": ["Hypertrophic cardiomyopathy"], "inheritance": "AD"},
    "TPM1": {"mim": "191010", "phenotypes": ["Hypertrophic cardiomyopathy"], "inheritance": "AD"},
    "MYL2": {"mim": "160781", "phenotypes": ["Hypertrophic cardiomyopathy"], "inheritance": "AD"},
    "MYL3": {"mim": "160790", "phenotypes": ["Hypertrophic cardiomyopathy"], "inheritance": "AD"},
    "ACTC1": {"mim": "102540", "phenotypes": ["Hypertrophic / dilated cardiomyopathy"], "inheritance": "AD"},
    "PRKAG2": {
        "mim": "602743",
        "phenotypes": ["Wolff-Parkinson-White / glycogen storage cardiomyopathy"],
        "inheritance": "AD",
    },
    "LMNA": {
        "mim": "150330",
        "phenotypes": ["Dilated cardiomyopathy", "Emery-Dreifuss muscular dystrophy"],
        "inheritance": "AD",
    },
    "KCNQ1": {"mim": "607542", "phenotypes": ["Long QT syndrome 1"], "inheritance": "AD"},
    "KCNH2": {"mim": "152427", "phenotypes": ["Long QT syndrome 2"], "inheritance": "AD"},
    "SCN5A": {"mim": "600163", "phenotypes": ["Long QT syndrome 3", "Brugada syndrome"], "inheritance": "AD"},
    "RYR2": {
        "mim": "180902",
        "phenotypes": ["Catecholaminergic polymorphic ventricular tachycardia"],
        "inheritance": "AD",
    },
    # --- Connective tissue / aortopathy -----------------------------------
    "FBN1": {"mim": "134797", "phenotypes": ["Marfan syndrome"], "inheritance": "AD"},
    "TGFBR1": {"mim": "190181", "phenotypes": ["Loeys-Dietz syndrome 1"], "inheritance": "AD"},
    "TGFBR2": {"mim": "190182", "phenotypes": ["Loeys-Dietz syndrome 2"], "inheritance": "AD"},
    "SMAD3": {"mim": "603109", "phenotypes": ["Loeys-Dietz syndrome 3"], "inheritance": "AD"},
    "ACTA2": {"mim": "102620", "phenotypes": ["Familial thoracic aortic aneurysm"], "inheritance": "AD"},
    "MYH11": {"mim": "160745", "phenotypes": ["Familial thoracic aortic aneurysm"], "inheritance": "AD"},
    "COL3A1": {"mim": "120180", "phenotypes": ["Ehlers-Danlos syndrome (vascular type)"], "inheritance": "AD"},
    # --- Familial hypercholesterolemia ------------------------------------
    "LDLR": {"mim": "606945", "phenotypes": ["Familial hypercholesterolemia"], "inheritance": "AD"},
    "APOB": {"mim": "107730", "phenotypes": ["Familial hypercholesterolemia (APOB)"], "inheritance": "AD"},
    "PCSK9": {"mim": "607786", "phenotypes": ["Familial hypercholesterolemia (PCSK9)"], "inheritance": "AD"},
    # --- RASopathies (Noonan-spectrum) ------------------------------------
    "PTPN11": {
        "mim": "176876",
        "phenotypes": ["Noonan syndrome 1", "LEOPARD syndrome", "Juvenile myelomonocytic leukemia"],
        "inheritance": "AD",
    },
    "SOS1": {"mim": "182530", "phenotypes": ["Noonan syndrome 4"], "inheritance": "AD"},
    "RAF1": {"mim": "164760", "phenotypes": ["Noonan syndrome 5", "LEOPARD syndrome"], "inheritance": "AD"},
    "KRAS": {
        "mim": "190070",
        "phenotypes": ["Noonan syndrome 3", "Cardio-facio-cutaneous syndrome"],
        "inheritance": "AD",
    },
    "BRAF": {
        "mim": "164757",
        "phenotypes": ["Cardio-facio-cutaneous syndrome", "LEOPARD syndrome"],
        "inheritance": "AD",
    },
    # --- Neurodevelopmental / chromatin (DDG2P-overlapping) ---------------
    "CHD8": {
        "mim": "610528",
        "phenotypes": ["Autism spectrum disorder", "Intellectual disability"],
        "inheritance": "AD",
    },
    "CHD7": {"mim": "608892", "phenotypes": ["CHARGE syndrome"], "inheritance": "AD"},
    "ARID1A": {"mim": "603024", "phenotypes": ["Coffin-Siris syndrome 2"], "inheritance": "AD"},
    "ARID1B": {"mim": "614556", "phenotypes": ["Coffin-Siris syndrome 1"], "inheritance": "AD"},
    "KMT2D": {"mim": "602113", "phenotypes": ["Kabuki syndrome 1"], "inheritance": "AD"},
    "KMT2A": {"mim": "159555", "phenotypes": ["Wiedemann-Steiner syndrome"], "inheritance": "AD"},
    "MECP2": {"mim": "300005", "phenotypes": ["Rett syndrome"], "inheritance": "XL"},
    "FMR1": {"mim": "309550", "phenotypes": ["Fragile X syndrome"], "inheritance": "XL"},
    "SCN1A": {
        "mim": "182389",
        "phenotypes": ["Dravet syndrome", "Generalized epilepsy with febrile seizures"],
        "inheritance": "AD",
    },
    "SCN2A": {
        "mim": "182390",
        "phenotypes": ["Developmental and epileptic encephalopathy", "Autism spectrum disorder"],
        "inheritance": "AD",
    },
    "KCNQ2": {
        "mim": "602235",
        "phenotypes": ["Benign familial neonatal seizures", "Epileptic encephalopathy"],
        "inheritance": "AD",
    },
    "STXBP1": {"mim": "602926", "phenotypes": ["Developmental and epileptic encephalopathy 4"], "inheritance": "AD"},
    "CASZ1": {"mim": "609895", "phenotypes": ["Cardiomyopathy / neurodevelopmental disorder"], "inheritance": "AD"},
    # --- Hematology / hemoglobinopathy ------------------------------------
    "HBB": {"mim": "141900", "phenotypes": ["Beta-thalassemia", "Sickle cell disease"], "inheritance": "AR"},
    "HBA1": {"mim": "141800", "phenotypes": ["Alpha-thalassemia"], "inheritance": "AR"},
    "HBA2": {"mim": "141850", "phenotypes": ["Alpha-thalassemia"], "inheritance": "AR"},
    "F8": {"mim": "300841", "phenotypes": ["Hemophilia A"], "inheritance": "XL"},
    "F9": {"mim": "300746", "phenotypes": ["Hemophilia B"], "inheritance": "XL"},
    # --- Renal / metabolic / pulmonary ------------------------------------
    "CFTR": {"mim": "602421", "phenotypes": ["Cystic fibrosis", "CBAVD"], "inheritance": "AR"},
    "PKD1": {"mim": "601313", "phenotypes": ["Polycystic kidney disease 1"], "inheritance": "AD"},
    "PKD2": {"mim": "173910", "phenotypes": ["Polycystic kidney disease 2"], "inheritance": "AD"},
    # --- Neurological / muscular ------------------------------------------
    "DMD": {"mim": "300377", "phenotypes": ["Duchenne / Becker muscular dystrophy"], "inheritance": "XL"},
    "SMN1": {"mim": "600354", "phenotypes": ["Spinal muscular atrophy"], "inheritance": "AR"},
    "HTT": {"mim": "613004", "phenotypes": ["Huntington disease"], "inheritance": "AD"},
    # --- Adult onset / complex --------------------------------------------
    "APOE": {
        "mim": "107741",
        "phenotypes": ["Alzheimer disease 2", "Hyperlipoproteinemia type III"],
        "inheritance": "AD/complex",
    },
    # --- Pharmacogenomics (legacy entries, retained for backward compat) --
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
