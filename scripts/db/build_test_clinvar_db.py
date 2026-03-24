"""Build a tiny ClinVar test DB with demo variants."""
import sqlite3
from pathlib import Path


def build_test_db(db_path: str = "data/db/clinvar_test.sqlite3"):
    Path(db_path).parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(db_path)

    conn.execute("""CREATE TABLE IF NOT EXISTS variants (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        chrom TEXT, pos INTEGER, ref TEXT, alt TEXT, rsid TEXT, gene TEXT,
        clinical_significance TEXT, review_status TEXT, phenotype_list TEXT,
        variation_id TEXT, allele_id TEXT, origin TEXT, assembly TEXT,
        last_evaluated TEXT, number_submitters INTEGER
    )""")
    conn.execute("""CREATE TABLE IF NOT EXISTS metadata (key TEXT PRIMARY KEY, value TEXT)""")

    variants = [
        ("chr17", 7675088, "C", "A", "rs28934578", "TP53", "Likely pathogenic", "reviewed by expert panel", "Li-Fraumeni syndrome", "182963", "182963", "germline", "GRCh38", "2024-09-06", 5),
        ("chr13", 32356550, "C", "T", "rs80358981", "BRCA2", "Uncertain significance", "criteria provided, single submitter", "Hereditary cancer", "801122", "801122", "germline", "GRCh38", "2023-01-01", 1),
        ("chr7", 117559590, "ATCT", "A", "rs113993960", "CFTR", "Pathogenic", "reviewed by expert panel", "Cystic fibrosis", "634837", "634837", "germline", "GRCh38", "2024-01-01", 10),
        ("chr10", 94781859, "G", "A", "rs4244285", "CYP2C19", "drug response", "criteria provided, multiple submitters, no conflicts", "CYP2C19 metabolism", "1703261", "1703261", "germline", "GRCh38", "2023-06-01", 8),
        ("chr6", 31464003, "T", "G", "rs2395029", "HLA-B", "Conflicting classifications of pathogenicity", "criteria provided, conflicting classifications", "Drug hypersensitivity", "14910", "14910", "germline", "GRCh38", "2022-01-01", 3),
        ("chr13", 48045719, "C", "T", "rs116855232", "NUDT15", "drug response", "criteria provided, multiple submitters, no conflicts", "Thiopurine response", "225201", "225201", "germline", "GRCh38", "2023-01-01", 4),
        ("chr2", 165342465, "G", "A", "rs794727152", "MUTYH", "Pathogenic", "criteria provided, single submitter", "MUTYH-associated polyposis", "1451854", "1451854", "germline", "GRCh38", "2023-06-01", 1),
        ("chr17", 72122759, "G", "A", "rs137853130", "PALB2", "Pathogenic", "reviewed by expert panel", "Hereditary breast cancer", "2517", "2517", "germline", "GRCh38", "2024-01-01", 6),
        ("chr12", 25225628, "C", "A", "rs121913527", "PTPN11", "Pathogenic/Likely pathogenic", "criteria provided, multiple submitters, no conflicts", "Noonan syndrome", "375963", "375963", "germline/somatic", "GRCh38", "2023-01-01", 4),
        ("chr19", 44908684, "T", "C", "rs429358", "APOE", "risk factor", "criteria provided, multiple submitters, no conflicts", "Alzheimer disease", "694742", "694742", "germline", "GRCh38", "2023-01-01", 3),
    ]

    conn.execute("DELETE FROM variants")
    conn.executemany("""INSERT INTO variants (chrom,pos,ref,alt,rsid,gene,clinical_significance,
        review_status,phenotype_list,variation_id,allele_id,origin,assembly,last_evaluated,number_submitters)
        VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""", variants)

    conn.execute("CREATE INDEX IF NOT EXISTS idx_chrom_pos ON variants(chrom, pos)")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_rsid ON variants(rsid)")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_chrom_pos_ref_alt ON variants(chrom, pos, ref, alt)")

    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('build_date', '2026-03-24T00:00:00')")
    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('source', 'test_data')")
    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('variant_count', '10')")
    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('assembly', 'GRCh38')")
    conn.execute("INSERT OR REPLACE INTO metadata VALUES ('clinvar_release', '2026-03')")

    conn.commit()
    conn.close()
    return db_path


if __name__ == "__main__":
    build_test_db()
