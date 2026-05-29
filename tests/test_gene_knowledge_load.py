"""gene_knowledge loader robustness (v2.7 ENRI-04).

_load_knowledge previously caught only FileNotFoundError / JSONDecodeError, so a
structurally-valid-but-wrong JSON (missing the top-level "genes" key, or entries
without a "gene" key) raised a KeyError that crashed the whole pipeline. It must
degrade to an empty mapping instead.
"""

import json

import scripts.common.gene_knowledge as gk


def test_load_knowledge_graceful_on_missing_genes_key(tmp_path, monkeypatch):
    bad = tmp_path / "bad.json"
    bad.write_text(json.dumps({"not_genes": []}))  # valid JSON, wrong schema
    monkeypatch.setattr(gk, "_KNOWLEDGE", None)
    monkeypatch.setattr(gk, "get", lambda *a, **k: str(bad))

    assert gk._load_knowledge() == {}  # no KeyError
    assert gk.get_gene_info("TP53") is None


def test_load_knowledge_skips_entries_without_gene_key(tmp_path, monkeypatch):
    f = tmp_path / "k.json"
    f.write_text(json.dumps({"genes": [{"gene": "TP53", "summary": "x"}, {"no_gene": True}]}))
    monkeypatch.setattr(gk, "_KNOWLEDGE", None)
    monkeypatch.setattr(gk, "get", lambda *a, **k: str(f))

    knowledge = gk._load_knowledge()
    assert "TP53" in knowledge
    assert len(knowledge) == 1  # the entry without a 'gene' key is skipped, not a crash


def test_load_knowledge_graceful_on_malformed_json(tmp_path, monkeypatch):
    bad = tmp_path / "bad.json"
    bad.write_text("{not valid json")
    monkeypatch.setattr(gk, "_KNOWLEDGE", None)
    monkeypatch.setattr(gk, "get", lambda *a, **k: str(bad))

    assert gk._load_knowledge() == {}
