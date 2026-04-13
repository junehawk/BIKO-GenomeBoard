"""Render Clinical Board opinion as HTML for report integration."""

from scripts.clinical_board.models import BoardOpinion


def render_board_opinion_html(opinion: BoardOpinion) -> str:
    """Render BoardOpinion as an HTML section for inclusion in the report."""
    if not opinion:
        return ""

    html_parts = []

    # Section header
    html_parts.append("""
    <div style="page-break-before:always;"></div>
    <div class="section-header">
      <span class="section-badge" style="background:#4338CA;">AI Clinical Board</span>
      <div class="section-rule" style="background:linear-gradient(90deg,#4338CA,#7C3AED);"></div>
    </div>

    <div style="background:#FEF3C7;border:1px solid #FCD34D;border-radius:8px;padding:12px 16px;margin-bottom:20px;font-size:12px;color:#92400E;">
      <strong>[AI-Generated]</strong> {disclaimer}
    </div>
    """.format(disclaimer=opinion.disclaimer))

    # Primary diagnosis
    html_parts.append(f"""
    <div style="background:#EEF2FF;border:1px solid #C7D2FE;border-radius:10px;padding:20px;margin-bottom:20px;">
      <div style="font-size:11px;font-weight:600;color:#6366F1;text-transform:uppercase;letter-spacing:0.5px;margin-bottom:6px;">Primary Diagnosis</div>
      <div style="font-size:18px;font-weight:700;color:#1E1B4B;margin-bottom:8px;">{opinion.primary_diagnosis or 'Not determined'}</div>
      <div style="font-size:13px;color:#4338CA;">{opinion.primary_diagnosis_evidence or ''}</div>
      <div style="margin-top:8px;">
        <span style="display:inline-block;background:{'#059669' if opinion.confidence == 'high' else '#D97706' if opinion.confidence == 'moderate' else '#DC2626'};color:#fff;border-radius:4px;padding:2px 10px;font-size:11px;font-weight:600;">
          Confidence: {opinion.confidence.upper()}
        </span>
        <span style="display:inline-block;background:#6366F1;color:#fff;border-radius:4px;padding:2px 10px;font-size:11px;font-weight:600;margin-left:6px;">
          Consensus: {opinion.agent_consensus}
        </span>
      </div>
    </div>
    """)

    # Key findings
    if opinion.key_findings:
        html_parts.append('<div style="margin-bottom:20px;">')
        html_parts.append('<div style="font-size:13px;font-weight:700;color:#1E293B;margin-bottom:8px;">Key Findings</div>')
        html_parts.append('<ul style="margin:0;padding-left:20px;font-size:13px;color:#334155;line-height:1.8;">')
        for f in opinion.key_findings:
            html_parts.append(f'<li>{f}</li>')
        html_parts.append('</ul></div>')

    # Differential diagnoses
    if opinion.differential_diagnoses:
        html_parts.append("""
        <div style="margin-bottom:20px;">
          <div style="font-size:13px;font-weight:700;color:#1E293B;margin-bottom:8px;">Differential Diagnoses</div>
          <table style="width:100%;border-collapse:collapse;font-size:12.5px;">
            <thead>
              <tr style="background:#F1F5F9;">
                <th style="text-align:left;padding:8px 12px;border-bottom:2px solid #E2E8F0;font-weight:600;">Diagnosis</th>
                <th style="text-align:center;padding:8px 12px;border-bottom:2px solid #E2E8F0;font-weight:600;width:100px;">Likelihood</th>
                <th style="text-align:left;padding:8px 12px;border-bottom:2px solid #E2E8F0;font-weight:600;">Evidence</th>
              </tr>
            </thead>
            <tbody>
        """)
        for dx in opinion.differential_diagnoses:
            likelihood = dx.get("likelihood", "unknown")
            color = "#059669" if likelihood == "high" else "#D97706" if likelihood == "moderate" else "#9CA3AF"
            html_parts.append(f"""
              <tr>
                <td style="padding:8px 12px;border-bottom:1px solid #E2E8F0;font-weight:600;">{dx.get('diagnosis', '')}</td>
                <td style="padding:8px 12px;border-bottom:1px solid #E2E8F0;text-align:center;">
                  <span style="color:{color};font-weight:600;">{likelihood}</span>
                </td>
                <td style="padding:8px 12px;border-bottom:1px solid #E2E8F0;color:#64748B;">{dx.get('evidence', '')}</td>
              </tr>
            """)
        html_parts.append('</tbody></table></div>')

    # Recommendations
    if opinion.recommendations:
        html_parts.append('<div style="margin-bottom:20px;">')
        html_parts.append('<div style="font-size:13px;font-weight:700;color:#1E293B;margin-bottom:8px;">Recommendations</div>')
        html_parts.append('<ol style="margin:0;padding-left:20px;font-size:13px;color:#334155;line-height:1.8;">')
        for r in opinion.recommendations:
            html_parts.append(f'<li>{r}</li>')
        html_parts.append('</ol></div>')

    # Follow-up
    if opinion.follow_up:
        html_parts.append('<div style="margin-bottom:20px;">')
        html_parts.append('<div style="font-size:13px;font-weight:700;color:#1E293B;margin-bottom:8px;">Follow-up</div>')
        html_parts.append('<ul style="margin:0;padding-left:20px;font-size:13px;color:#64748B;line-height:1.8;">')
        for f in opinion.follow_up:
            html_parts.append(f'<li>{f}</li>')
        html_parts.append('</ul></div>')

    # Dissenting opinions
    if opinion.dissenting_opinions:
        html_parts.append("""
        <div style="background:#FFF7ED;border:1px solid #FED7AA;border-radius:8px;padding:12px 16px;margin-bottom:20px;">
          <div style="font-size:12px;font-weight:700;color:#C2410C;margin-bottom:6px;">Dissenting Opinions</div>
          <ul style="margin:0;padding-left:18px;font-size:12.5px;color:#9A3412;line-height:1.7;">
        """)
        for d in opinion.dissenting_opinions:
            html_parts.append(f'<li>{d}</li>')
        html_parts.append('</ul></div>')

    # Individual agent summaries (collapsible)
    if opinion.agent_opinions:
        html_parts.append("""
        <div style="margin-bottom:20px;">
          <div style="font-size:13px;font-weight:700;color:#1E293B;margin-bottom:10px;">Domain Specialist Opinions</div>
        """)
        for agent_op in opinion.agent_opinions:
            confidence_color = "#059669" if agent_op.confidence == "high" else "#D97706" if agent_op.confidence == "moderate" else "#DC2626"
            findings_html = ""
            for f in agent_op.findings[:3]:  # Show top 3 findings
                finding_text = f.get("finding", f) if isinstance(f, dict) else str(f)
                findings_html += f'<li>{finding_text}</li>'

            html_parts.append(f"""
            <details style="margin-bottom:8px;border:1px solid #E2E8F0;border-radius:6px;overflow:hidden;">
              <summary style="padding:10px 14px;background:#F8FAFC;cursor:pointer;font-size:13px;font-weight:600;color:#1E293B;">
                {agent_op.agent_name}
                <span style="float:right;color:{confidence_color};font-size:11px;font-weight:600;">{agent_op.confidence.upper()}</span>
              </summary>
              <div style="padding:12px 14px;font-size:12.5px;color:#475569;">
                <ul style="margin:0 0 8px;padding-left:18px;line-height:1.7;">{findings_html}</ul>
                {'<div style="font-size:11px;color:#94A3B8;margin-top:4px;">Refs: ' + ', '.join(agent_op.references[:3]) + '</div>' if agent_op.references else ''}
              </div>
            </details>
            """)
        html_parts.append('</div>')

    return "\n".join(html_parts)
