#!/usr/bin/env python3
"""Build local gnomAD SQLite database from sites VCF files."""

import gzip
import sqlite3
import logging
import sys
from pathlib import Path
from datetime import datetime

logger = logging.getLogger(__name__)

DEFAULT_DB_PATH = "data/db/gnomad.sqlite3"


def build_db(vcf_paths: list, db_path: str = DEFAULT_DB_PATH, version: str = "4.1") -> str:
    """Build SQLite DB from gnomAD sites VCF file(s).

    gnomAD VCF INFO fields of interest:
    - AF: global allele frequency
    - AF_eas: East Asian allele frequency
    - AF_afr, AF_amr, AF_nfe, AF_sas, etc.
    - AN: total allele number
    - AC: total allele count
    """
    Path(db_path).parent.mkdir(parents=True, exist_ok=True)

    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA synchronous=NORMAL")

    conn.execute("""
        CREATE TABLE IF NOT EXISTS variants (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            chrom TEXT NOT NULL,
            pos INTEGER NOT NULL,
            ref TEXT NOT NULL,
            alt TEXT NOT NULL,
            rsid TEXT,
            af_global REAL,
            af_eas REAL,
            af_afr REAL,
            af_amr REAL,
            af_nfe REAL,
            af_sas REAL,
            an INTEGER,
            ac INTEGER,
            filter_status TEXT
        )
    """)

    conn.execute("""
        CREATE TABLE IF NOT EXISTS metadata (
            key TEXT PRIMARY KEY,
            value TEXT
        )
    """)

    total_count = 0

    for vcf_path in vcf_paths:
        logger.info(f"Processing {vcf_path}...")
        count = 0
        batch = []

        open_func = gzip.open if str(vcf_path).endswith('.gz') else open
        mode = 'rt' if str(vcf_path).endswith('.gz') else 'r'

        with open_func(vcf_path, mode) as f:
            for line in f:
                if line.startswith('#'):
                    continue

                fields = line.strip().split('\t')
                if len(fields) < 8:
                    continue

                chrom = fields[0]
                if not chrom.startswith('chr'):
                    chrom = f"chr{chrom}"

                try:
                    pos = int(fields[1])
                except ValueError:
                    continue

                rsid = fields[2] if fields[2] != '.' else None
                ref = fields[3]
                alt_field = fields[4]
                filter_status = fields[6]
                info = fields[7]

                # Parse multiple ALT alleles
                for alt_idx, alt in enumerate(alt_field.split(',')):
                    # Parse INFO field for frequencies
                    af_global = None
                    af_eas = None
                    af_afr = None
                    af_amr = None
                    af_nfe = None
                    af_sas = None
                    an = None
                    ac = None

                    for item in info.split(';'):
                        if '=' not in item:
                            continue
                        key, val = item.split('=', 1)

                        # Handle multi-allelic: pick the value for this alt allele
                        vals = val.split(',')
                        v = vals[alt_idx] if alt_idx < len(vals) else vals[0]

                        try:
                            if key == 'AF':
                                af_global = float(v)
                            elif key == 'AF_eas':
                                af_eas = float(v)
                            elif key == 'AF_afr':
                                af_afr = float(v)
                            elif key == 'AF_amr':
                                af_amr = float(v)
                            elif key == 'AF_nfe':
                                af_nfe = float(v)
                            elif key == 'AF_sas':
                                af_sas = float(v)
                            elif key == 'AN':
                                an = int(v)
                            elif key == 'AC':
                                ac = int(vals[alt_idx] if alt_idx < len(vals) else vals[0])
                        except (ValueError, IndexError):
                            pass

                    batch.append((
                        chrom, pos, ref, alt, rsid,
                        af_global, af_eas, af_afr, af_amr, af_nfe, af_sas,
                        an, ac, filter_status
                    ))
                    count += 1

                if len(batch) >= 10000:
                    conn.executemany("""
                        INSERT INTO variants (chrom, pos, ref, alt, rsid,
                            af_global, af_eas, af_afr, af_amr, af_nfe, af_sas,
                            an, ac, filter_status)
                        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    """, batch)
                    batch = []
                    if count % 500000 == 0:
                        logger.info(f"  {count:,} variants...")

        if batch:
            conn.executemany("""
                INSERT INTO variants (chrom, pos, ref, alt, rsid,
                    af_global, af_eas, af_afr, af_amr, af_nfe, af_sas,
                    an, ac, filter_status)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, batch)

        total_count += count
        logger.info(f"  {vcf_path}: {count:,} variants")

    # Create indexes
    logger.info("Creating indexes...")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_chrom_pos ON variants(chrom, pos)")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_chrom_pos_ref_alt ON variants(chrom, pos, ref, alt)")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_rsid ON variants(rsid)")

    # Metadata
    now = datetime.utcnow().isoformat()
    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('build_date', ?)", (now,))
    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('gnomad_version', ?)", (version,))
    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('variant_count', ?)", (str(total_count),))
    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('assembly', 'GRCh38')")

    conn.commit()
    conn.close()

    logger.info(f"gnomAD DB built: {total_count:,} variants → {db_path}")
    return db_path


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
    if len(sys.argv) < 2:
        print("Usage: python build_gnomad_db.py <vcf1.gz> [vcf2.gz ...] [--version 4.1] [--output path]")
        sys.exit(1)

    vcf_paths = [p for p in sys.argv[1:] if not p.startswith('--')]
    version = "4.1"
    db_path = DEFAULT_DB_PATH
    for i, arg in enumerate(sys.argv):
        if arg == '--version' and i + 1 < len(sys.argv):
            version = sys.argv[i + 1]
        if arg == '--output' and i + 1 < len(sys.argv):
            db_path = sys.argv[i + 1]

    build_db(vcf_paths, db_path, version)
