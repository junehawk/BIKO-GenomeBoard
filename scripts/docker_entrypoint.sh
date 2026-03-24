#!/bin/bash
# GenomeBoard Docker convenience commands

case "$1" in
  single)
    shift
    python scripts/orchestrate.py "$@"
    ;;
  batch)
    shift
    python scripts/orchestrate.py --batch "$@"
    ;;
  build-clinvar-db)
    shift
    python scripts/db/build_clinvar_db.py "$@"
    ;;
  build-gnomad-db)
    shift
    python scripts/db/build_gnomad_db.py "$@"
    ;;
  test)
    pip install -r requirements-dev.txt
    python -m pytest tests/ -v
    ;;
  *)
    python scripts/orchestrate.py "$@"
    ;;
esac
