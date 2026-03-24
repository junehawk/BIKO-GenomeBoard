"""Query gnomAD using tabix-indexed VCF files (pysam).

This module queries downloaded gnomAD .vcf.bgz + .tbi files directly,
without any SQLite conversion. Files are expected at:
  data/db/gnomad_vcf/gnomad.exomes.v4.1.sites.chr{N}.vcf.bgz

Each chromosome has its own file + .tbi index.
"""

import logging
import os
import re
import threading
from pathlib import Path
from typing import Dict, List, Optional

from scripts.common.config import get
from scripts.common.models import Variant

logger = logging.getLogger(__name__)

_tabix_files = {}  # chrom -> pysam.TabixFile
_lock = threading.Lock()

def _get_vcf_dir() -> str:
    return get("paths.gnomad_vcf_dir", "data/db/gnomad_vcf")

def _find_vcf_for_chrom(chrom: str) -> Optional[str]:
    """Find the gnomAD VCF file for a given chromosome."""
    vcf_dir = _get_vcf_dir()
    # Normalize: chr17 -> 17, chrX -> X
    chrom_num = chrom.replace("chr", "")

    # Try common gnomAD file patterns
    patterns = [
        f"gnomad.exomes.v4.1.sites.chr{chrom_num}.vcf.bgz",
        f"gnomad.exomes.v4.1.sites.{chrom_num}.vcf.bgz",
        f"gnomad.genomes.v4.1.sites.chr{chrom_num}.vcf.bgz",
        f"gnomad.exomes.*.chr{chrom_num}.vcf.bgz",
    ]

    for pattern in patterns:
        # Direct match first
        path = os.path.join(vcf_dir, pattern)
        if os.path.exists(path):
            return path

    # Glob fallback
    vcf_dir_path = Path(vcf_dir)
    if vcf_dir_path.exists():
        for f in vcf_dir_path.glob(f"*chr{chrom_num}*.vcf.bgz"):
            if f.with_suffix(".bgz.tbi").exists() or Path(str(f) + ".tbi").exists():
                return str(f)

    return None

def _get_tabix(chrom: str):
    """Get or open a pysam.TabixFile for a chromosome."""
    with _lock:
        if chrom in _tabix_files:
            return _tabix_files[chrom]

    vcf_path = _find_vcf_for_chrom(chrom)
    if not vcf_path:
        return None

    tbi_path = vcf_path + ".tbi"
    if not os.path.exists(tbi_path):
        logger.warning(f"Tabix index not found: {tbi_path}")
        return None

    try:
        import pysam
        tbx = pysam.TabixFile(vcf_path)
        with _lock:
            _tabix_files[chrom] = tbx
        return tbx
    except Exception as e:
        logger.warning(f"Failed to open tabix file {vcf_path}: {e}")
        return None

def _parse_info_field(info: str) -> Dict[str, str]:
    """Parse VCF INFO field into dict."""
    result = {}
    for item in info.split(";"):
        if "=" in item:
            key, val = item.split("=", 1)
            result[key] = val
        else:
            result[item] = "true"
    return result

def _extract_af(info_dict: Dict[str, str], alt_idx: int = 0) -> Dict:
    """Extract allele frequencies from gnomAD INFO fields."""
    def _get_af(key):
        val = info_dict.get(key)
        if val is None:
            return None
        try:
            parts = val.split(",")
            return float(parts[alt_idx] if alt_idx < len(parts) else parts[0])
        except (ValueError, IndexError):
            return None

    return {
        "af_global": _get_af("AF"),
        "af_eas": _get_af("AF_eas"),
        "af_afr": _get_af("AF_afr"),
        "af_amr": _get_af("AF_amr"),
        "af_nfe": _get_af("AF_nfe"),
        "af_sas": _get_af("AF_sas"),
    }

def query_tabix_gnomad(variant: Variant) -> Dict:
    """Query gnomAD via tabix for a single variant.

    Returns same structure as query_gnomad() API and query_local_gnomad() SQLite.
    """
    chrom_num = variant.chrom.replace("chr", "")
    tbx = _get_tabix(variant.chrom)

    if tbx is None:
        return {"gnomad_all": None, "gnomad_eas": None, "api_available": False}

    try:
        # Tabix query: fetch records overlapping this position
        # gnomAD VCFs use "chr1" format for GRCh38
        query_chrom = variant.chrom if variant.chrom.startswith("chr") else f"chr{variant.chrom}"

        # Try with and without chr prefix
        records = None
        for qc in [query_chrom, chrom_num]:
            try:
                records = list(tbx.fetch(qc, variant.pos - 1, variant.pos))
                if records:
                    break
            except ValueError:
                continue

        if not records:
            return {"gnomad_all": None, "gnomad_eas": None, "api_available": True}

        # Find exact match (pos + ref + alt)
        for record in records:
            fields = record.split("\t")
            if len(fields) < 8:
                continue

            rec_pos = int(fields[1])
            rec_ref = fields[3]
            rec_alts = fields[4].split(",")
            rec_info = fields[7]

            if rec_pos != variant.pos or rec_ref != variant.ref:
                continue

            # Check each ALT allele
            for alt_idx, alt in enumerate(rec_alts):
                if alt == variant.alt:
                    info_dict = _parse_info_field(rec_info)
                    afs = _extract_af(info_dict, alt_idx)

                    return {
                        "gnomad_all": afs["af_global"],
                        "gnomad_eas": afs["af_eas"],
                        "gnomad_afr": afs.get("af_afr"),
                        "gnomad_amr": afs.get("af_amr"),
                        "gnomad_nfe": afs.get("af_nfe"),
                        "gnomad_sas": afs.get("af_sas"),
                        "api_available": True,
                    }

        # Position found but no exact alt match
        return {"gnomad_all": None, "gnomad_eas": None, "api_available": True}

    except Exception as e:
        logger.warning(f"Tabix query failed for {variant.variant_id}: {e}")
        return {"gnomad_all": None, "gnomad_eas": None, "api_available": False}

def get_available_chromosomes() -> List[str]:
    """List which chromosomes have gnomAD VCF files available."""
    vcf_dir = Path(_get_vcf_dir())
    if not vcf_dir.exists():
        return []
    chroms = []
    for f in sorted(vcf_dir.glob("*.vcf.bgz")):
        m = re.search(r'chr(\w+)', f.name)
        if m:
            chroms.append(f"chr{m.group(1)}")
    return chroms

def get_db_version() -> Dict:
    """Get gnomAD tabix metadata."""
    vcf_dir = Path(_get_vcf_dir())
    if not vcf_dir.exists():
        return {"source": "not_available"}

    files = list(vcf_dir.glob("*.vcf.bgz"))
    if not files:
        return {"source": "not_available"}

    # Infer version from filename
    version = "unknown"
    for f in files:
        m = re.search(r'v(\d+\.\d+)', f.name)
        if m:
            version = m.group(1)
            break

    # Get modification date of first file
    mod_date = files[0].stat().st_mtime
    from datetime import datetime

    return {
        "source": "tabix",
        "gnomad_version": version,
        "build_date": datetime.fromtimestamp(mod_date).strftime("%Y-%m-%d"),
        "assembly": "GRCh38",
        "chromosomes": len(files),
        "vcf_dir": str(vcf_dir),
    }

def close():
    """Close all open tabix files."""
    with _lock:
        for tbx in _tabix_files.values():
            try:
                tbx.close()
            except Exception:
                pass
        _tabix_files.clear()
