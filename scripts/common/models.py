# scripts/common/models.py
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, List, Optional
import re

STRENGTH_MAP = {
    "PVS": "very_strong",
    "PS": "strong",
    "PM": "moderate",
    "PP": "supporting",
    "BA": "stand_alone",
    "BS": "strong",
    "BP": "supporting",
}


@dataclass
class Variant:
    chrom: str
    pos: int
    ref: str
    alt: str
    gene: Optional[str] = None
    rsid: Optional[str] = None  # rs number from VCF ID column
    # VEP/SnpEff annotation fields (from pre-annotated VCF)
    hgvsc: Optional[str] = None  # c. notation (e.g., c.524G>T)
    hgvsp: Optional[str] = None  # p. notation (e.g., p.Arg175Leu)
    consequence: Optional[str] = None  # e.g., missense_variant, frameshift_variant
    transcript: Optional[str] = None  # e.g., NM_000546.6 or ENST00000269305
    impact: Optional[str] = None  # HIGH, MODERATE, LOW, MODIFIER
    sift: Optional[str] = None  # e.g., deleterious(0.01)
    polyphen: Optional[str] = None  # e.g., probably_damaging(0.998)
    in_silico: Optional[Dict] = None  # REVEL, CADD, AlphaMissense, SpliceAI scores

    @classmethod
    def from_string(cls, s: str) -> Variant:
        """Parse 'chr17:7577120 G>A' or '17:7577120 G>A'"""
        m = re.match(r"(chr)?(\w+):(\d+)\s+([ACGT]+)>([ACGT]+)", s.strip())
        if not m:
            raise ValueError(f"Cannot parse variant: {s}")
        chrom = f"chr{m.group(2)}" if not m.group(1) else f"chr{m.group(2)}"
        return cls(chrom=chrom, pos=int(m.group(3)), ref=m.group(4), alt=m.group(5))

    @property
    def variant_id(self) -> str:
        return f"{self.chrom}:{self.pos}:{self.ref}>{self.alt}"


@dataclass
class AcmgEvidence:
    code: str
    source: str
    description: str

    @property
    def strength(self) -> str:
        base = re.match(r"([A-Z]+)", self.code)
        if base:
            prefix = base.group(1)
            if "_Supporting" in self.code:
                return "supporting"
            return STRENGTH_MAP.get(prefix, "unknown")
        return "unknown"

    @property
    def direction(self) -> str:
        if self.code.startswith(("PVS", "PS", "PM", "PP")):
            return "pathogenic"
        return "benign"


@dataclass
class FrequencyData:
    krgdb: Optional[float] = None
    gnomad_eas: Optional[float] = None
    gnomad_all: Optional[float] = None
    korea4k: Optional[float] = None
    nard2: Optional[float] = None

    @property
    def korean_max(self) -> Optional[float]:
        """Maximum frequency across all Korean population sources."""
        korean_freqs = [f for f in [self.krgdb, self.korea4k, self.nard2] if f is not None]
        return max(korean_freqs) if korean_freqs else None

    def korean_vs_global_ratio(self) -> Optional[float]:
        kr = self.korean_max
        if kr is not None and self.gnomad_all and self.gnomad_all > 0:
            return kr / self.gnomad_all
        return None


@dataclass
class PgxResult:
    gene: str
    star_allele: str
    phenotype: str
    cpic_level: str
    korean_prevalence: float
    western_prevalence: float
    clinical_impact: str = ""
    cpic_recommendation: str = ""

    @property
    def korean_flag(self) -> bool:
        if self.western_prevalence > 0:
            return (self.korean_prevalence / self.western_prevalence) >= 2.0
        return self.korean_prevalence > 0


@dataclass
class StructuralVariant:
    annotsv_id: str
    chrom: str
    start: int
    end: int
    length: int
    sv_type: str              # DEL, DUP, INV, BND, INS
    sample_id: str
    acmg_class: int           # 1-5
    ranking_score: float
    cytoband: str
    gene_name: str            # "ERBB2" or "TBX1;COMT;HIRA"
    gene_count: int
    # Gene detail (populated from split rows)
    gene_details: List[Dict] = field(default_factory=list)
    # Pathogenic evidence
    p_gain_phen: str = ""
    p_gain_hpo: str = ""
    p_gain_source: str = ""
    p_loss_phen: str = ""
    p_loss_hpo: str = ""
    p_loss_source: str = ""
    # Benign evidence
    b_gain_af_max: Optional[float] = None
    b_loss_af_max: Optional[float] = None
    # OMIM
    omim_morbid: bool = False

    @property
    def is_pathogenic(self) -> bool:
        return self.acmg_class in (4, 5)

    @property
    def is_vus(self) -> bool:
        return self.acmg_class == 3

    @property
    def is_benign(self) -> bool:
        return self.acmg_class in (1, 2)

    @property
    def acmg_label(self) -> str:
        labels = {5: "Pathogenic", 4: "Likely Pathogenic", 3: "VUS", 2: "Likely Benign", 1: "Benign"}
        return labels.get(self.acmg_class, "Unknown")

    @property
    def size_display(self) -> str:
        if abs(self.length) >= 1_000_000:
            return f"{abs(self.length) / 1_000_000:.1f} Mb"
        elif abs(self.length) >= 1000:
            return f"{abs(self.length) / 1000:.1f} kb"
        return f"{abs(self.length)} bp"

    @property
    def phenotypes(self) -> str:
        if self.sv_type == "DUP" and self.p_gain_phen:
            return self.p_gain_phen
        if self.sv_type == "DEL" and self.p_loss_phen:
            return self.p_loss_phen
        return self.p_gain_phen or self.p_loss_phen or ""

    @property
    def evidence_source(self) -> str:
        if self.sv_type == "DUP" and self.p_gain_source:
            return self.p_gain_source
        if self.sv_type == "DEL" and self.p_loss_source:
            return self.p_loss_source
        return self.p_gain_source or self.p_loss_source or ""

    def is_dosage_sensitive(self, mode: str = "cancer") -> bool:
        """Check if this Class 3 VUS is dosage-sensitive enough to display."""
        if self.acmg_class != 3:
            return False
        hi_thresh = 1 if mode == "rare-disease" else 2
        pli_thresh = 0.8 if mode == "rare-disease" else 0.9
        for gd in self.gene_details:
            hi = gd.get("hi") or 0
            ts = gd.get("ts") or 0
            pli = gd.get("pli") or 0.0
            if self.sv_type == "DEL" and hi >= hi_thresh:
                return True
            if self.sv_type == "DUP" and ts >= hi_thresh:
                return True
            if pli >= pli_thresh:
                return True
        # Large VUS with multiple OMIM genes
        if abs(self.length) > 1_000_000 and self.gene_count >= 3 and self.omim_morbid:
            return True
        return False
