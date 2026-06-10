"""HTML-escaping of LLM-derived Clinical Board fields (T3-01 / v2.7 D3).

The board HTML fragment is emitted into the report via Jinja {{ ...|safe }},
which bypasses autoescape. Every LLM/curated free-text value interpolated into
that fragment must therefore be html.escape()'d at the point of interpolation,
or a model emission containing markup becomes live HTML in the report.
"""

from scripts.clinical_board.models import BoardOpinion, CancerBoardOpinion
from scripts.clinical_board.render import render_board_opinion_html

PAYLOAD = "<script>alert(1)</script>"
ESCAPED = "&lt;script&gt;"


def test_cancer_opinion_escapes_llm_fields():
    opinion = CancerBoardOpinion(
        therapeutic_headline=PAYLOAD,
        therapeutic_implications=PAYLOAD,
        therapeutic_evidence=PAYLOAD,
        treatment_options=[
            {"drug": "<img src=x onerror=alert(1)>", "evidence_level": "A", "resistance_notes": PAYLOAD}
        ],
        actionable_findings=[PAYLOAD],
        clinical_actions=[PAYLOAD],
        immunotherapy_eligibility=PAYLOAD,
        monitoring_plan=[PAYLOAD],
        dissenting_opinions=[PAYLOAD],
    )
    out = render_board_opinion_html(opinion, language="en")

    # The literal LLM markup must never survive into the rendered fragment.
    assert PAYLOAD not in out
    assert "<img src=x onerror=alert(1)>" not in out
    # It must instead appear escaped.
    assert ESCAPED in out
    assert "&lt;img src=x" in out


def test_rare_disease_opinion_escapes_llm_fields():
    opinion = BoardOpinion(
        primary_diagnosis=PAYLOAD,
        primary_diagnosis_evidence=PAYLOAD,
        key_findings=[PAYLOAD],
        recommendations=[PAYLOAD],
        differential_diagnoses=[{"diagnosis": PAYLOAD, "likelihood": "high"}],
        follow_up=[PAYLOAD],
        dissenting_opinions=[PAYLOAD],
    )
    out = render_board_opinion_html(opinion, language="en")

    assert PAYLOAD not in out
    assert ESCAPED in out


def test_agent_opinion_findings_are_escaped():
    from scripts.clinical_board.models import AgentOpinion

    opinion = BoardOpinion(
        primary_diagnosis="benign",
        agent_opinions=[
            AgentOpinion(
                agent_name="<script>Pathologist</script>",
                domain="variant_pathology",
                findings=[{"finding": PAYLOAD}],
                references=[PAYLOAD],
            )
        ],
    )
    out = render_board_opinion_html(opinion, language="en")

    assert PAYLOAD not in out
    assert "<script>Pathologist</script>" not in out
    assert ESCAPED in out
