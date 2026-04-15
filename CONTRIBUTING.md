# Contributing to BIKO GenomeBoard

Thank you for your interest in contributing. BIKO GenomeBoard is a
**research reference tool** that produces genomic interpretation
documents for independent review by a researcher or clinician. It is
**not** a clinical decision instrument and is not cleared for direct
patient care — please keep this framing in mind when proposing
features or filing issues.

## Prerequisites

- **Python** 3.10, 3.11, or 3.12 (CI runs all three).
- Ability to run `setup_databases.sh` to build the local reference
  databases (ClinVar, gnomAD, CIViC, HPO, OMIM, ClinGen, Korean pop
  sources, etc.). Some tests are DB-state-independent; integration
  tests are not.
- **Optional**: a local Ollama install if you plan to work on the
  Clinical Board LLM path. The pipeline runs without it — the
  deterministic curator + template renderer fallback carry the
  evidence path.

BIKO uses only local databases and Ollama — there are **no external
API keys** required to run the core pipeline.

## Development setup

```bash
git clone git@github.com:junehawk/BIKO-GenomeBoard.git
cd BIKO-GenomeBoard
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
pip install -r requirements-dev.txt
```

Then build the local databases:

```bash
bash setup_databases.sh
```

## Pre-push checklist — run these BEFORE `git push`

CI enforces **ruff `0.15.10`** for lint and format, plus `pytest` across
Python 3.10 / 3.11 / 3.12. The #1 cause of CI red on this repo is a
missed `ruff format` pass on an edited file: unit tests pass locally,
`ruff check` passes locally, but `ruff format --check` fails on CI
because the editor-saved formatting differs from ruff's canonical
output. Run the full four-step sequence locally before every push:

```bash
# 1. Lint — autofix what ruff can, review the rest
ruff check --fix scripts/ tests/

# 2. Format — must match CI's `ruff format --check` exactly
ruff format scripts/ tests/

# 3. Verify both gates are green (what CI actually runs)
ruff check scripts/ tests/           # must say "All checks passed!"
ruff format --check scripts/ tests/  # must say "N files already formatted"

# 4. Full test suite (CI runs without the integration marker)
python -m pytest tests/ -q
```

The full checklist (with rationale, historical notes, and edge cases
around consequence-form normalization and DB-state-independent tests)
lives in [`CLAUDE.md`](./CLAUDE.md) under **Pre-push checklist — run
these BEFORE `git push`**. Please read it once before your first PR.

## Running tests

```bash
python -m pytest tests/ -q              # what CI runs (non-integration)
python -m pytest tests/ -x              # stop on first failure
python -m pytest tests/test_foo.py -v   # single file, verbose
```

Tests that depend on optional local DBs (`civic.sqlite3`,
`hpo.sqlite3`, etc.) must guard on **actual query results**, not on
`os.path.exists` — CI runners do not provision these and an empty stub
can trip existence checks while queries return nothing.

## Pull request workflow

1. **Branch off `main`.** One focused change per PR — no drive-by
   refactors bundled into a feature change.
2. **Conventional Commits** for commit messages. This repo already
   follows this style, so please match it:
   - `feat: ...` — new capability
   - `fix: ...` — bug fix
   - `docs: ...` — documentation only
   - `test: ...` — test-only change
   - `chore: ...` — tooling, config, hygiene
   - `ci: ...` — CI configuration
3. Run the **Pre-push checklist** (above) before pushing.
4. Open the PR against `main`. In the description, explain the *why*,
   link any related issue, and note any follow-up work you chose to
   defer.
5. Do **not** push directly to `main`. Do **not** use `--no-verify`
   or skip hooks.

## Privacy — synthetic data only

**Do not** add real patient data, PHI, real VCFs from actual samples,
or any clinically-identifying information to this repository — not in
code, tests, issues, PR descriptions, commit messages, or
documentation. Synthetic test data only.

If you discover PHI or patient-identifying content anywhere in the
repository (including git history), follow the private-reporting
process in [`SECURITY.md`](./SECURITY.md). Do not discuss the finding
in a public issue or PR.

Several paths are excluded via `.gitignore` specifically to keep
working material and sample inputs out of the tracked tree:
`_workspace/`, `docs/superpowers/plans/`, `data/sample_vcf/codegen-*`,
`docs/API_KEYS.md`. Please do not reintroduce them.

## Roadmap and deferred work

The current deferred work and v2.3+ roadmap is documented in the
release notes for each version — see [`CHANGELOG.md`](./CHANGELOG.md)
and the v2.2.0 / v2.2.1 GitHub release bodies. The short version:
Board Chair LLM schema work, PM5 wiring in `classify.py`, MSI /
fusion / mutational signatures, PMID references on agent opinions,
and patient-header metadata persistence are all on the v2.3+ list.

Architectural do-not-touch zones are enumerated under **Things to
never do** in [`CLAUDE.md`](./CLAUDE.md) — notably, the classification
and tiering engine must stay deterministic Python and never sit
behind an LLM call.

## Questions

For general questions, open a GitHub Discussion or a feature-request
issue. For anything involving PHI, security, or a vulnerability,
follow [`SECURITY.md`](./SECURITY.md) instead.
