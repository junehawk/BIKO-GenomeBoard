# Security Policy

## Research-use-only disclaimer

BIKO GenomeBoard produces a **research reference document**. It is **not**
a clinical decision instrument, not a diagnostic device, and is not
cleared for direct patient care. Security incidents in this repository
therefore have research-data implications rather than direct
patient-care implications. Any downstream clinical use is the
responsibility of the reviewing researcher or clinician, who makes
their own independent judgement.

## Supported versions

| Version       | Supported                              |
| ------------- | -------------------------------------- |
| `2.2.x`       | ✅ Current — security fixes accepted   |
| `< 2.2.0`     | ❌ Unsupported — please upgrade        |

The rollback checkpoint `pre-v2.2-phaseA` is preserved as a git tag for
continuity only and is **not** a supported baseline.

## Reporting a vulnerability

Please report suspected vulnerabilities **privately**, not via a public
GitHub issue. Two channels are accepted:

1. **GitHub Security Advisories** (preferred):
   <https://github.com/junehawk/BIKO-GenomeBoard/security/advisories/new>
2. **Email**: `[REDACTED — fill in before publishing]`

Include the affected version, a minimal reproduction (synthetic data
only — see below), and the impact you observed. Please allow
reasonable time for triage before any public disclosure.

## PHI / patient-data incidents

If you discover patient-identifying information, real genomic data, or
any other PHI in a commit, issue, pull request, or release artifact:

- **Do not** discuss the finding in any public venue, including GitHub
  issues, PR comments, discussions, or external channels.
- **Do not** fork, mirror, quote, or redistribute the affected content.
- Contact the maintainer privately through one of the channels above
  and clearly label the report as a **PHI incident**.

Repository hygiene is enforced through `.gitignore` rules that exclude
internal working material and sample inputs that may carry sensitive
content, including:

- `_workspace/`
- `docs/superpowers/plans/`
- `data/sample_vcf/codegen-*`
- `docs/API_KEYS.md`

If you observe any of these paths reappearing in a commit, flag it as a
PHI incident.

## Response window

This project is a researcher-maintained repository without a dedicated
security team. Reports are handled on a **best-effort** basis and there
is **no formal SLA**. We will acknowledge reports as promptly as is
practical and keep you informed as remediation progresses.

## Safe harbor

Good-faith, responsible disclosure is welcomed. If you follow this
policy — private reporting, synthetic reproductions, no public
discussion before remediation — we will not pursue or support any
action against you for the research involved.
