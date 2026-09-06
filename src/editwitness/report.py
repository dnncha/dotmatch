"""Self-contained, script-free HTML. All reference names, descriptions and sequences are escaped."""
from __future__ import annotations

from html import escape

from .models import Analysis, Assay

MAX_REPORT_WITNESSES = 50
MAX_SEQUENCE_PREVIEW_BASES = 2400
MAX_REPORT_PRODUCTS = 8
MAX_REPORT_SEQUENCE_BASES = 250_000

CSS = """
:root{color-scheme:light;--ink:#142832;--muted:#506671;--line:#d8e2e5;--paper:#f5f8f9;
--accent:#126d72;--warn:#855014}*{box-sizing:border-box}body{margin:0;background:var(--paper);
color:var(--ink);overflow-wrap:anywhere;font:16px/1.6 -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif}
header{background:#112b36;color:#fff;padding:42px max(24px,calc((100vw - 1060px)/2)) 34px}
.brand{font-size:14px;letter-spacing:.12em;font-weight:750;text-transform:uppercase;color:#9ce2d8}
h1{font-size:clamp(28px,4vw,42px);line-height:1.15;letter-spacing:-.035em;margin:16px 0}
header p{max-width:820px;color:#d4e3e9;margin-bottom:0}.tag{display:inline-block;font-size:12px;
font-weight:700;letter-spacing:.06em;text-transform:uppercase;border:1px solid #557783;
border-radius:4px;padding:3px 8px;margin-left:12px;color:#fff}main{max-width:1108px;margin:auto;
padding:26px 24px 48px}.banner{border-left:4px solid var(--warn);padding:14px 20px;
background:#fff6e8;border-radius:0 7px 7px 0;margin:0 0 24px}.grid{display:grid;
grid-template-columns:repeat(3,1fr);gap:16px}.metric{padding:20px 22px;background:white;
border:1px solid var(--line);border-radius:10px}.metric strong{display:block;font-size:36px;
line-height:1.2;letter-spacing:-.04em}.metric span{color:var(--muted);font-size:14px}
section{background:#fff;border:1px solid var(--line);border-radius:10px;padding:24px;margin-top:24px}
h2{font-size:23px;line-height:1.3;letter-spacing:-.02em;margin:0 0 12px}h3{font-size:17px;
margin:0 0 8px}p{margin:10px 0}.muted,.small{color:var(--muted)}.small{font-size:13px}
code,pre{font-family:ui-monospace,SFMono-Regular,Consolas,monospace;font-size:13px}
code{overflow-wrap:anywhere}pre{background:#f3f6f7;padding:14px;border-radius:6px;
white-space:pre-wrap;overflow-wrap:anywhere;max-height:260px;overflow:auto}table{border-collapse:collapse;
width:100%;font-size:14px}caption{text-align:left;font-size:13px;color:var(--muted);padding-bottom:12px}
th,td{padding:12px 10px;text-align:left;vertical-align:top;border-bottom:1px solid var(--line)}
th{font-size:12px;text-transform:uppercase;letter-spacing:.04em;background:#f7f9fa}
.scroll{overflow-x:auto}.scroll:focus-visible{outline:3px solid var(--accent);outline-offset:3px}
table{min-width:620px}table code{white-space:nowrap}details{border-top:1px solid var(--line);padding:14px 0}
summary{cursor:pointer;font-weight:650;overflow-wrap:anywhere}summary:focus-visible{outline:3px solid var(--accent)}
.witness{border-left:3px solid #ae722c;padding:0 0 0 16px;margin:22px 0}
.pill{display:inline-block;border-radius:4px;padding:2px 7px;font-size:12px;font-weight:650;
background:#edf2f4;color:#344e5c}.pill.warn{background:#fff0d7;color:#805019}
svg{display:block;width:100%;height:auto;min-width:940px}ul,ol{padding-left:23px}li{margin:8px 0}
footer{padding-top:24px;color:var(--muted);font-size:12px;overflow-wrap:anywhere}
@media(max-width:650px){header{padding:28px 20px}.grid{grid-template-columns:1fr}.metric{display:flex;
align-items:center;gap:16px}.metric strong{font-size:30px}main{padding:20px 16px}section{padding:18px}
th,td{padding:10px 8px}.tag{margin-left:6px;font-size:10px}}
@media print{body{background:#fff;font-size:11pt}header{background:#fff;color:#142832;padding:0 0 16px}
header p,.brand{color:#142832}.tag{color:#142832}main{padding:0}.grid{display:flex}.metric{flex:1}
section{break-inside:avoid;page-break-inside:avoid;border-color:#aaa}.scroll{overflow:visible}
pre{max-height:none}details{display:block}footer{font-size:8pt}}
"""


def e(value: object) -> str:
    return escape(str(value), quote=True)


def _assay_map(assays: tuple[Assay, ...], length: int) -> str:
    height = 55 + 52 * len(assays)
    scale = 730 / length
    parts = [f'<svg viewBox="0 0 940 {height}" role="img" aria-label="Annotated primer positions on the local reference">',
             '<text x="175" y="18" fill="#506671" font-size="12">Local reference coordinates · 0-based, half-open</text>']
    for index, assay in enumerate(assays):
        y = 45 + index * 52
        parts += [f'<text x="0" y="{y+5}" fill="#142832" font-size="13"><title>{e(assay.id)}</title>{e(assay.id[:20])}{"…" if len(assay.id) > 20 else ""}</text>',
                  f'<line x1="175" y1="{y}" x2="905" y2="{y}" stroke="#d4dfe3" stroke-width="3"/>']
        for interval, label in ((assay.left_primer, "L"), (assay.right_primer, "R")):
            x = 175 + interval.start * scale
            width = max(2, (interval.end - interval.start) * scale)
            parts += [f'<rect x="{x:.3f}" y="{y-7}" width="{width:.3f}" height="14" rx="2" fill="#126d72"/>',
                      f'<text x="{x + width if label == "R" else x:.3f}" text-anchor="{"end" if label == "R" else "start"}" y="{y+24}" font-size="12" fill="#506671">{label} [{interval.start},{interval.end})</text>']
    parts.append('</svg>')
    return "".join(parts)


def render_report(result: Analysis) -> str:
    """Render model output only; no scientific computations are performed here."""
    n = len(result.witnesses)
    unresolved = len(result.plan.unresolved_hypotheses)
    title = "Your assay leaves alternatives unresolved" if n else "No alternative matched within the declared model"
    intro = (
        f"{n} declared alternative {'hypotheses produce' if n != 1 else 'hypothesis produces'} the same modeled observations "
        f"as {e(result.expected_hypothesis)}. This demonstrates ambiguity, not an observed editing defect."
        if n else "Only the supplied hypotheses were tested. This result does not establish complete validation "
                  "or exclude biological outcomes that were not modeled."
    )
    rows = []
    for h in result.hypotheses:
        expected = h.hypothesis_id == result.expected_hypothesis
        identical = h.same_local_genomic_state_as_expected and not expected
        status = "Expected hypothesis" if expected else (
            "Same local DNA state; excluded" if identical else (
                "Same observations" if h.equivalent_to_expected else "Different observations"
            )
        )
        css = "pill warn" if h.equivalent_to_expected and not expected and not identical else "pill"
        rows.append(f'<tr><td><code>{e(h.hypothesis_id)}</code></td>'
                    f'<td>{e(" + ".join(h.alleles))}</td><td><span class="{css}">{status}</span></td>'
                    f'<td>{e(", ".join(h.distinguishing_existing_assays)) or "—"}</td></tr>')
    witnesses = []
    preview_remaining = MAX_REPORT_SEQUENCE_BASES

    def read_preview(reads: tuple[str, ...]) -> str:
        nonlocal preview_remaining
        parts = []
        for index, read in enumerate(reads):
            limit = min(MAX_SEQUENCE_PREVIEW_BASES, preview_remaining)
            if len(read) and not limit:
                parts.append(f"read {index + 1}: [report preview budget reached; full sequence in JSON]")
                continue
            shown = read[:limit]
            preview_remaining -= len(shown)
            suffix = (f" … [preview: first {limit} of {len(read)} bases; full sequence in JSON]"
                      if len(read) > limit else "")
            parts.append(f"read {index + 1}: {shown or '(empty insert)'}{suffix}")
        return "\n".join(parts)

    by_pair = {(o.allele_id, o.assay_id): o for o in result.allele_observations}
    for w in result.witnesses[:MAX_REPORT_WITNESSES]:
        evidence = []
        ids = set(w.expected_alleles + w.alternative_alleles)
        for allele in result.allele_evidence:
            if allele.allele_id not in ids:
                continue
            edits = []
            for edit in allele.edits:
                action = "delete" if not edit.sequence else ("insert" if edit.start == edit.end else "replace")
                replacement = (" → " + edit.sequence[:120]
                               + (f" … [{len(edit.sequence)} bases; full definition in JSON]"
                                  if len(edit.sequence) > 120 else "")) if edit.sequence else ""
                edits.append(f"{action} reference[{edit.start}:{edit.end}){replacement}")
            evidence.append(
                f'<h3>Allele definition: {e(allele.allele_id)}</h3><p>{e(allele.description)}</p>'
                f'<p class="small">Final local sequence: {allele.sequence_length:,} bp. '
                f'SHA-256: <code>{e(allele.sequence_sha256)}</code></p>'
                f'<pre>{e(chr(10).join(edits) or "Unchanged local reference sequence")}</pre>'
            )
        relevant_assays = result.assays + tuple(
            a for a in result.candidates if a.id in w.resolving_candidate_assays
        )
        for assay in relevant_assays:
            for allele_id in sorted(ids):
                observation = by_pair.get((allele_id, assay.id))
                if observation is None:
                    continue
                role = "existing" if assay in result.assays else "separating candidate"
                evidence.append(f'<h3>{e(assay.id)} / {e(allele_id)} · {role}</h3>'
                    f'<p class="small">{e(observation.status)}. {e(observation.reason)}</p>')
                if observation.reads:
                    evidence.append(f'<pre>{e(read_preview(observation.reads))}</pre>')
                elif observation.products:
                    for product in observation.products[:MAX_REPORT_PRODUCTS]:
                        evidence.append(f'<p class="small">{e(product.orientation)} product: '
                            f'final-allele sites [{product.plus_left_site.start}:{product.plus_left_site.end}) / '
                            f'[{product.plus_right_site.start}:{product.plus_right_site.end}); '
                            f'{product.product_length} bp.</p>'
                            f'<pre>{e(read_preview(product.reads))}</pre>')
                    if len(observation.products) > MAX_REPORT_PRODUCTS:
                        evidence.append(f'<p class="small">First {MAX_REPORT_PRODUCTS} of '
                            f'{len(observation.products)} products shown; every product is in the JSON.</p>')
                else:
                    evidence.append('<p class="small">No modeled sequence signal.</p>')
        witnesses.append(
            f'<article class="witness"><h3>{e(w.hypothesis_id)}</h3><p><code>{e(" + ".join(w.expected_alleles))}</code>'
            f' versus <code>{e(" + ".join(w.alternative_alleles))}</code></p><p>{e(w.explanation)}</p>'
            f'<p class="small">Separating candidates: {e(", ".join(w.resolving_candidate_assays)) or "None supplied"}.</p>'
            f'<details><summary>Inspect allele definitions and sequence evidence</summary>{"".join(evidence)}</details></article>'
        )
    report_limit = (f"<p class=\"small\">Showing the first {MAX_REPORT_WITNESSES} of {n} counterexamples. "
                    "The JSON result contains every declared counterexample.</p>"
                    if n > MAX_REPORT_WITNESSES else "")
    generation_note = ""
    if result.generation is not None:
        g = result.generation
        generation_note = (f"<p class=\"small\">Declared generation: {g.valid_deletions:,} grid deletions; "
                           f"{g.added_hypotheses:,} hypotheses added; {g.deduplicated_states:,} duplicate "
                           "local sequence states omitted. These are design-grid counts, not outcome frequencies.</p>")
    choices = e(", ".join(result.plan.selected_assays)) or "None"
    plan_note = (
        f'<p><strong>{unresolved} alternative(s) remain unresolved by every supplied candidate.</strong> '
        'Repeating the same sequence-presence measurement does not add dosage evidence. '
        'Consider an independently justified orthogonal measurement; this alpha does not simulate it.</p>'
        if unresolved else '<p>Candidate coverage applies only to the declared alternatives and idealized response model.</p>'
    )
    notices = "".join(f'<li><code>{e(x.code)}</code> — {e(x.message)}'
                      f'{" [" + e(", ".join(x.related_ids)) + "]" if x.related_ids else ""}</li>'
                      for x in result.notices)
    assumptions = "".join(f'<li>{e(item)}</li>' for item in result.assumptions)
    assay_rows = "".join(
        f'<tr><td><code>{e(a.id)}</code></td><td>{e(a.readout)}'
        f'{"; " + str(a.read_bases) + " insert bases/end" if a.read_bases is not None else ""}</td>'
        f'<td>{a.min_product_bp}–{a.max_product_bp if a.max_product_bp is not None else "unbounded"} bp</td>'
        f'<td>{a.cost_units}</td></tr>' for a in result.assays + result.candidates
    )
    return f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<meta http-equiv="Content-Security-Policy" content="default-src 'none'; style-src 'unsafe-inline'; base-uri 'none'; form-action 'none'">
<meta name="referrer" content="no-referrer"><title>EditWitness · {e(result.reference.name)}</title>
<style>{CSS}</style></head><body>
<header><div class="brand">EditWitness <span class="tag">Research alpha</span></div>
<h1>{title}</h1><p>{intro}</p><p class="small" style="color:#d4e3e9">Observation model: <code>{e(result.model_version)}</code></p></header><main>
<div class="banner"><strong>Model evidence, not a safety certificate.</strong> Software-tested; not empirically validated.
{' This is a synthetic demonstration.' if result.reference.synthetic else ''}</div>
<div class="grid"><div class="metric"><strong>{n}</strong><span>equivalent alternative hypotheses</span></div>
<div class="metric"><strong>{len(result.plan.resolved_hypotheses)}</strong><span>separable with supplied candidates</span></div>
<div class="metric"><strong>{unresolved}</strong><span>beyond the supplied candidates</span></div></div>
<section><h2>What was compared</h2><p class="small">Reference: {e(result.reference.name)} · {result.reference.length:,} bp.
 Sequence presence only; no inference from read fractions or allele dosage.</p>{generation_note}<div class="scroll" tabindex="0" role="region" aria-label="Scrollable evidence table or assay map">
<table><caption>Existing assays only. “Different observations” is conditional on the response model.</caption>
<thead><tr><th>Hypothesis</th><th>Alleles</th><th>Comparison</th><th>Separating existing assays</th></tr></thead>
<tbody>{''.join(rows)}</tbody></table></div></section>
<section><h2>Concrete counterexamples</h2>{report_limit}{''.join(witnesses) or '<p>No matching alternative was declared. That is not proof of completeness.</p>'}</section>
<section><h2>Which additional assays help?</h2><p>Selected panel: <strong>{choices}</strong> · {result.plan.cost_units} declared cost units.</p>
<p class="small">Method: {e(result.plan.algorithm)}. Optimality: {e(result.plan.optimality)}.</p>
<p>{e(result.plan.note)}</p>{plan_note}</section>
<section><h2>Assay geometry and readout</h2><div class="scroll" tabindex="0" role="region" aria-label="Scrollable evidence table or assay map">{_assay_map(result.assays + result.candidates, result.reference.length)}</div>
<div class="scroll" tabindex="0" role="region" aria-label="Scrollable evidence table or assay map"><table><caption>Cost units are user supplied; product-size bounds are declared inclusion rules.</caption>
<thead><tr><th>Assay</th><th>Observed sequence</th><th>Product bounds</th><th>Cost</th></tr></thead>
<tbody>{assay_rows}</tbody></table></div></section>
<section><h2>Important qualifications</h2><ul>{notices}</ul></section>
<section><h2>The model behind this report</h2><ol>{assumptions}</ol></section>
<footer>EditWitness {e(result.package_version)} · {e(result.model_version)}<br>
Normalized manifest SHA-256: <code>{e(result.manifest_sha256)}</code><br>
Result SHA-256: <code>{e(result.result_sha256)}</code><br>
This report is local, self-contained and script-free. It may contain sensitive genomic sequence; share deliberately.
 A checksum checks content integrity, not the truth of a scientific claim.</footer>
</main></body></html>'''
