"""Build a tiny gnomAD test DB with demo variants."""
import sqlite3
from pathlib import Path


def build_test_db(db_path: str = "data/db/gnomad_test.sqlite3"):
    Path(db_path).parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(db_path)

    conn.execute("""CREATE TABLE IF NOT EXISTS variants (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        chrom TEXT, pos INTEGER, ref TEXT, alt TEXT, rsid TEXT,
        af_global REAL, af_eas REAL, af_afr REAL, af_amr REAL, af_nfe REAL, af_sas REAL,
        an INTEGER, ac INTEGER, filter_status TEXT
    )""")
    conn.execute("""CREATE TABLE IF NOT EXISTS metadata (key TEXT PRIMARY KEY, value TEXT)""")

    variants = [
        ("chr17", 7675088, "C", "A", "rs28934578", 0.000004, 0.000008, 0.000001, 0.0, 0.000005, 0.0, 1520000, 6, "PASS"),
        ("chr13", 32356550, "C", "T", "rs80358981", 0.0000131, 0.0000182, 0.0000095, 0.0, 0.0000146, 0.0, 1520000, 20, "PASS"),
        ("chr7", 117559590, "ATCT", "A", "rs113993960", 0.00788, 0.0001, 0.001, 0.003, 0.017, 0.003, 1520000, 11978, "PASS"),
        ("chr10", 94781859, "G", "A", "rs4244285", 0.168, 0.30, 0.15, 0.12, 0.15, 0.17, 1520000, 255360, "PASS"),
        ("chr6", 31464003, "T", "G", "rs2395029", 0.0247, 0.00853, 0.003, 0.01, 0.05, 0.02, 1520000, 37544, "PASS"),
        ("chr13", 48045719, "C", "T", "rs116855232", 0.01257, 0.101, 0.001, 0.002, 0.001, 0.003, 1520000, 19107, "PASS"),
    ]

    conn.execute("DELETE FROM variants")
    conn.executemany("""INSERT INTO variants (chrom,pos,ref,alt,rsid,af_global,af_eas,af_afr,af_amr,af_nfe,af_sas,an,ac,filter_status)
        VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)""", variants)

    conn.execute("CREATE INDEX IF NOT EXISTS idx_chrom_pos ON variants(chrom, pos)")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_rsid ON variants(rsid)")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_chrom_pos_ref_alt ON variants(chrom, pos, ref, alt)")

    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('build_date', '2026-03-24T00:00:00')")
    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('gnomad_version', '4.1')")
    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('variant_count', '6')")
    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('assembly', 'GRCh38')")

    conn.commit()
    conn.close()
    return db_path


if __name__ == "__main__":
    build_test_db()
