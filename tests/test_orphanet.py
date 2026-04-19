# tests/test_orphanet.py
import sqlite3
from pathlib import Path

import pytest


def _create_sample_orphanet_xml(path):
    """Minimal Orphanet Product 9 XML for testing."""
    Path(path).write_text("""<?xml version="1.0" encoding="UTF-8"?>
<JDBOR date="2025-12-01" version="1.0">
  <DisorderList count="2">
    <Disorder id="1">
      <OrphaCode>166024</OrphaCode>
      <Name lang="en">Cystic fibrosis</Name>
      <DisorderType><Name lang="en">Disease</Name></DisorderType>
      <DisorderGeneAssociationList count="1">
        <DisorderGeneAssociation>
          <Gene><Symbol>CFTR</Symbol><Name lang="en">CF transmembrane conductance regulator</Name></Gene>
        </DisorderGeneAssociation>
      </DisorderGeneAssociationList>
      <PrevalenceList count="1">
        <Prevalence>
          <PrevalenceType><Name lang="en">Point prevalence</Name></PrevalenceType>
          <PrevalenceQualification><Name lang="en">Value and class</Name></PrevalenceQualification>
          <PrevalenceClass><Name lang="en">1-5 / 10 000</Name></PrevalenceClass>
          <ValMoy>3.0</ValMoy>
          <PrevalenceGeographic><Name lang="en">Europe</Name></PrevalenceGeographic>
        </Prevalence>
      </PrevalenceList>
    </Disorder>
    <Disorder id="2">
      <OrphaCode>558</OrphaCode>
      <Name lang="en">Marfan syndrome</Name>
      <DisorderType><Name lang="en">Disease</Name></DisorderType>
      <DisorderGeneAssociationList count="1">
        <DisorderGeneAssociation>
          <Gene><Symbol>FBN1</Symbol><Name lang="en">fibrillin 1</Name></Gene>
        </DisorderGeneAssociation>
      </DisorderGeneAssociationList>
      <PrevalenceList count="1">
        <Prevalence>
          <PrevalenceType><Name lang="en">Point prevalence</Name></PrevalenceType>
          <PrevalenceQualification><Name lang="en">Value and class</Name></PrevalenceQualification>
          <PrevalenceClass><Name lang="en">1-9 / 100 000</Name></PrevalenceClass>
          <ValMoy>1.5</ValMoy>
          <PrevalenceGeographic><Name lang="en">Worldwide</Name></PrevalenceGeographic>
        </Prevalence>
      </PrevalenceList>
    </Disorder>
  </DisorderList>
</JDBOR>""")


def test_build_orphanet_db(tmp_path):
    from scripts.db.build_orphanet_db import build_db

    xml_path = str(tmp_path / "orphanet.xml")
    db_path = str(tmp_path / "orphanet.sqlite3")
    _create_sample_orphanet_xml(xml_path)
    result = build_db(xml_path, db_path)
    assert result == db_path
    conn = sqlite3.connect(db_path)
    rows = conn.execute("SELECT * FROM prevalence WHERE gene_symbol = 'CFTR'").fetchall()
    assert len(rows) >= 1
    conn.close()


def test_query_prevalence_by_gene(tmp_orphanet_db):
    from scripts.db.query_orphanet import get_prevalence_by_gene

    result = get_prevalence_by_gene("CFTR", tmp_orphanet_db)
    assert len(result) >= 1
    assert result[0]["disease_name"] == "Cystic fibrosis"
    assert "1-5 / 10 000" in result[0]["prevalence_class"]


def test_get_prevalence_text(tmp_orphanet_db):
    from scripts.db.query_orphanet import get_prevalence_text

    text = get_prevalence_text("CFTR", tmp_orphanet_db)
    assert "Cystic fibrosis" in text
    assert "1-5 / 10 000" in text


def test_query_unknown_gene(tmp_orphanet_db):
    from scripts.db.query_orphanet import get_prevalence_by_gene

    result = get_prevalence_by_gene("FAKEGENE", tmp_orphanet_db)
    assert result == []


@pytest.fixture
def tmp_orphanet_db(tmp_path):
    from scripts.db.build_orphanet_db import build_db

    xml_path = str(tmp_path / "orphanet.xml")
    db_path = str(tmp_path / "orphanet.sqlite3")
    _create_sample_orphanet_xml(xml_path)
    build_db(xml_path, db_path)
    return db_path
