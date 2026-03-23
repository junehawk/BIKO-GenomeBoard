# scripts/common/models.py
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional
import re

STRENGTH_MAP = {
    "PVS": "very_strong",
    "PS": "strong", "PM": "moderate", "PP": "supporting",
    "BA": "stand_alone", "BS": "strong", "BP": "supporting",
}

@dataclass
class Variant:
    chrom: str
    pos: int
    ref: str
    alt: str
    gene: Optional[str] = None
    rsid: Optional[str] = None  # rs number from VCF ID column

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

    def korean_vs_global_ratio(self) -> Optional[float]:
        if self.krgdb is not None and self.gnomad_all and self.gnomad_all > 0:
            return self.krgdb / self.gnomad_all
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
