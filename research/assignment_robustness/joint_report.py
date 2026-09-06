#!/usr/bin/env python3
"""Validate complete AR003 outputs and analyse gene-annotation representation.

No phenotype p-values, gene discovery or unconditional confidence intervals.
The old and new implementations' complete read budgets are retained separately.
"""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import html
import json
import math
from pathlib import Path
import sqlite3
import sys
import tempfile
import unittest

POLICIES = ('exact', 'radius_k1', 'best_k1')
ACCESSIONS = ('ERR376998', 'ERR376999', 'SRR8297997')
STATES = ('unique', 'ambiguous', 'none', 'invalid')
AUDIT_RUN = 34032669590
AUDIT_COMMIT = '68a406b41805b808e6428927c2b8ccd2a2e55f9f'
PROTOCOL_COMMIT = 'de0002a37259cf27000015f43dd335343cd01234'


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n', encoding='utf-8')


def write_tsv(path, header, rows):
    with Path(path).open('w', newline='', encoding='utf-8') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(header)
        w.writerows(rows)


def read_tsv(path):
    with Path(path).open(newline='', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f, delimiter='\t')
        if not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError('Missing or repeated columns: '+str(path))
        rows = list(reader)
    if any(None in r or any(v is None for v in r.values()) for r in rows):
        raise ValueError('Nonrectangular table: '+str(path))
    return rows


def natural(text):
    if not isinstance(text, str) or not text.isascii() or not text.isdecimal():
        raise ValueError('Expected a nonnegative integer, not '+repr(text))
    return int(text)


def manifests(root):
    found = {}
    for p in sorted(Path(root).rglob('completion.json')):
        obj = json.loads(p.read_text())
        accession = obj.get('accession')
        if accession not in ACCESSIONS:
            continue
        if accession in found or obj.get('completion') != 'complete':
            raise ValueError('Duplicate or incomplete evidence for '+accession)
        for name, expected in obj['files'].items():
            item = p.parent / name
            if not item.resolve().is_relative_to(p.parent.resolve()) or not item.is_file() or digest(item) != expected:
                raise ValueError('Evidence checksum/path failed: '+str(item))
        found[accession] = (p.parent, obj, digest(p))
    if set(found) != set(ACCESSIONS):
        raise ValueError('Exactly three completed archives are required')
    return found


def keyed(rows, column):
    result = {}
    for row in rows:
        key = row[column]
        if not key or key in result:
            raise ValueError('Missing or duplicate identity: '+repr(key))
        result[key] = row
    return result


def recount(joint, old):
    """Independently reconstruct published counts/bounds from every class."""
    folder, manifest, _ = joint
    old_folder, previous, _ = old
    if (manifest['input']['sha256'] != previous['input']['sha256'] or
        manifest['library_sha256'] != previous['library_sha256'] or
        manifest['records'] != int(previous['input']['read_count']) or
        manifest['offsets'] != previous['selected_offsets']):
        raise AssertionError('Old/new input or search domain differs')
    rows = read_tsv(old_folder / 'all-guide-counts.tsv')
    ids = [r['id'] for r in rows]
    genes = [r['gene'] for r in rows]
    keyed(rows, 'id')
    guide_output = keyed(read_tsv(folder / 'guide-counts.tsv'), 'id')
    gene_output = keyed(read_tsv(folder / 'gene-counts-and-bounds.tsv'), 'gene')
    if set(guide_output) != set(ids) or set(gene_output) != set(genes):
        raise AssertionError('Reference identities were dropped or added')
    guide = [Counter() for _ in POLICIES]
    lower = [Counter() for _ in POLICIES]
    upper = [Counter() for _ in POLICIES]
    states = [[Counter(), Counter()] for _ in POLICIES]
    empties = Counter()
    seen = set()
    classes = read_tsv(folder / 'candidate-classes.tsv')
    for r in classes:
        p = POLICIES.index(r['policy'])
        text = r['target_indices']
        if (p, text) in seen:
            raise AssertionError('Repeated candidate class')
        seen.add((p, text))
        n = natural(r['records'])
        if n <= 0:
            raise ValueError('Empty-count candidate class')
        targets = tuple(natural(x) for x in text.split(',')) if text else ()
        if tuple(sorted(set(targets))) != targets or any(i >= len(rows) for i in targets):
            raise ValueError('Malformed candidate identity set')
        if not targets:
            empties[p] += n
            continue
        group = {genes[i] for i in targets}
        states[p][0]['unique' if len(targets) == 1 else 'ambiguous'] += n
        states[p][1]['unique' if len(group) == 1 else 'ambiguous'] += n
        if len(targets) == 1:
            guide[p][targets[0]] += n
        for g in group:
            upper[p][g] += n
            if len(group) == 1:
                lower[p][g] += n
    for p, policy in enumerate(POLICIES):
        for j, resolution in enumerate(('guide', 'gene')):
            q = manifest['qc'][policy][resolution]
            if any(states[p][j][s] != q[s] for s in ('unique', 'ambiguous')):
                raise AssertionError('Class/state reconstruction failed')
            if empties[p] != q['none'] + q['invalid'] or sum(q.values()) != manifest['records']:
                raise AssertionError('Unmatched/invalid budget failed')
        if any(natural(guide_output[ident][policy]) != guide[p][i] or guide_output[ident]['gene'] != genes[i] for i, ident in enumerate(ids)):
            raise AssertionError('Class/guide count reconstruction failed')
        if any(natural(gene_output[g][policy+'_lower']) != lower[p][g] or natural(gene_output[g][policy+'_upper']) != upper[p][g] for g in gene_output):
            raise AssertionError('Class/gene bounds reconstruction failed')
        if any(lower[p][g] > upper[p][g] for g in gene_output):
            raise AssertionError('Inverted count range')
        if sum(lower[p].values()) != manifest['qc'][policy]['gene']['unique']:
            raise AssertionError('Sum of lower counts does not equal gene-resolved budget')
    controls = json.loads((folder / 'known-origin-controls.json').read_text())
    if controls['records'] != len(rows):
        raise AssertionError('Control/reference population mismatch')
    for policy in POLICIES:
        q = controls['policies'][policy]
        n = controls['records']
        for resolution in ('guide', 'gene'):
            if sum(q.get(resolution+'_'+s, 0) for s in ('unique', 'ambiguous', 'none')) + q.get('invalid', 0) != n:
                raise AssertionError('Known-origin control budget failed')
            if q.get(resolution+'_correct', 0)+q.get(resolution+'_incorrect', 0) != q.get(resolution+'_unique', 0):
                raise AssertionError('Control correctness labels fail to partition unique calls')
    return rows, gene_output, controls, len(classes)


def legacy_groups(rows):
    """Python and SQL independently aggregate every original annotation."""
    result = defaultdict(lambda: dict(exact=0, radius_k1=0, best_k1=0, supported_guides=0))
    con = sqlite3.connect(':memory:')
    con.execute('CREATE TABLE counts(id TEXT PRIMARY KEY, gene TEXT, exact INTEGER, radius INTEGER, best INTEGER, supported INTEGER)')
    inserts = []
    for r in rows:
        numbers = [natural(r[p]) for p in POLICIES]
        support = int(numbers[0] >= 10)
        for p, n in zip(POLICIES, numbers):
            result[r['gene']][p] += n
        result[r['gene']]['supported_guides'] += support
        inserts.append((r['id'], r['gene'], *numbers, support))
    con.executemany('INSERT INTO counts VALUES(?,?,?,?,?,?)', inserts)
    sql = {g: (a, b, c, s) for g, a, b, c, s in con.execute('SELECT gene,SUM(exact),SUM(radius),SUM(best),SUM(supported) FROM counts GROUP BY gene')}
    con.close()
    if set(sql) != set(result) or any(sql[g] != tuple(result[g][k] for k in (*POLICIES, 'supported_guides')) for g in sql):
        raise AssertionError('Independent SQL/Python annotation aggregation differs')
    return dict(result)


def log_ratio(treatment, control, n_control, n_treatment, pseudo):
    if min(treatment, control) < 0 or min(n_control, n_treatment, pseudo) <= 0:
        raise ValueError('Invalid count, exposure or pseudocount')
    answer = math.log2((treatment+pseudo)/(control+pseudo) * n_control/n_treatment)
    check = math.log2(treatment+pseudo)-math.log2(control+pseudo)+math.log2(n_control)-math.log2(n_treatment)
    if not math.isfinite(answer) or abs(answer-check) > 2e-12:
        raise AssertionError('Independent ratio/log-difference arithmetic failed')
    return answer


def analyse_pair(old_rows, joint_gene, manifests_by_id):
    control = legacy_groups(old_rows['ERR376998'])
    treatment = legacy_groups(old_rows['ERR376999'])
    if set(control) != set(treatment):
        raise AssertionError('Control/treatment annotation population mismatch')
    nc = manifests_by_id['ERR376998'][1]['records']
    nt = manifests_by_id['ERR376999'][1]['records']
    gc = joint_gene['ERR376998']
    gt = joint_gene['ERR376999']
    if set(gc) != set(control) or set(gt) != set(control):
        raise AssertionError('Joint/fixed annotation population mismatch')
    detailed, rows, arithmetic_checks = [], [], 0
    summaries = []
    for threshold in (20, 50, 100):
        eligible = sorted(g for g in control if control[g]['exact'] >= threshold and control[g]['supported_guides'] >= 2)
        for pseudo in (0.1, 0.5, 1.0):
            primary = threshold == 50 and pseudo == 0.5
            for policy in POLICIES:
                shifts, sign_reversals, intervals_excluding_zero, deltas, widths = 0, 0, 0, [], []
                for g in sorted(control):
                    c, t = control[g]['best_k1'], treatment[g]['best_k1']
                    baseline = log_ratio(t, c, nc, nt, pseudo)
                    cl, cu = natural(gc[g][policy+'_lower']), natural(gc[g][policy+'_upper'])
                    tl, tu = natural(gt[g][policy+'_lower']), natural(gt[g][policy+'_upper'])
                    point = log_ratio(tl, cl, nc, nt, pseudo)
                    low = log_ratio(tl, cu, nc, nt, pseudo)
                    high = log_ratio(tu, cl, nc, nt, pseudo)
                    arithmetic_checks += 4
                    if not low-1e-12 <= point <= high+1e-12:
                        raise AssertionError('Point not within its conditional effect range')
                    delta = point-baseline
                    included = control[g]['exact'] >= threshold and control[g]['supported_guides'] >= 2
                    if primary:
                        detailed.append((g, policy, int(included), control[g]['exact'], control[g]['supported_guides'], c, t, cl, tl, cu, tu, baseline, point, delta, low, high))
                    if not included:
                        continue
                    deltas.append(abs(delta))
                    widths.append(high-low)
                    shifts += abs(delta) >= 0.5
                    sign_reversals += baseline*point < 0 and min(abs(baseline), abs(point)) >= 0.5
                    intervals_excluding_zero += low > 0 or high < 0
                if not eligible:
                    raise ValueError('No eligible genes; do not invent a percentage')
                deltas.sort()
                widths.sort()
                summary = {'policy': policy, 'baseline_exact_minimum': threshold, 'pseudocount': pseudo,
                           'primary': primary, 'eligible_annotations': len(eligible), 'shift_ge_half_log2': shifts,
                           'shift_fraction': shifts/len(eligible), 'strong_sign_reversals': sign_reversals,
                           'median_absolute_delta': deltas[(len(deltas)-1)//2], 'max_absolute_delta': max(deltas),
                           'median_conditional_range_width': widths[(len(widths)-1)//2],
                           'conditional_range_excludes_zero': intervals_excluding_zero,
                           'normalization': 'input_read_totals', 'reference': 'fixed_window_best_k1',
                           'biological_FDR_or_hit_calling': False}
                summaries.append(summary)
    return detailed, summaries, arithmetic_checks


def number(n):
    return format(n, ',') if isinstance(n, int) else format(n, '.6g')


def markdown_table(headers, rows):
    out = ['| '+' | '.join(headers)+' |', '|'+ '|'.join('---' for _ in headers)+'|']
    for row in rows:
        out.append('| '+' | '.join(str(v).replace('|', '\\|') for v in row)+' |')
    return '\n'.join(out)


def build(args):
    output = Path(args.out)
    if output.exists():
        raise FileExistsError(output)
    output.mkdir(parents=True, exist_ok=False)
    joint = manifests(args.joint)
    old = manifests(args.replay)
    old_rows, new_genes, controls, class_counts = {}, {}, {}, {}
    integrity = []
    for accession in ACCESSIONS:
        old_rows[accession], new_genes[accession], controls[accession], class_counts[accession] = recount(joint[accession], old[accession])
        for kind, values in (('joint', joint), ('prior_replay', old)):
            folder, manifest, checksum = values[accession]
            integrity.append({'accession': accession, 'kind': kind, 'completion_sha256': checksum, 'verified_files': len(manifest['files'])})
    if sum(joint[a][1]['records'] for a in ACCESSIONS) != 30215791:
        raise AssertionError('Locked full-archive population differs')
    detail, effects, arithmetic = analyse_pair(old_rows, new_genes, joint)
    qc_rows = []
    control_rows = []
    summary_samples = {}
    for accession in ACCESSIONS:
        obj = joint[accession][1]
        summary_samples[accession] = {k: obj[k] for k in ('records', 'target_count', 'gene_annotations', 'offsets', 'qc', 'additional_gene_resolved_records', 'validation', 'input', 'library_sha256', 'environment')}
        for policy in POLICIES:
            for resolution in ('guide', 'gene', 'position'):
                qc_rows.append((accession, policy, resolution, *(obj['qc'][policy][resolution][s] for s in STATES)))
            c = controls[accession]['policies'][policy]
            control_rows.append((accession, policy, controls[accession]['records'], c.get('origin_in_candidate_set', 0),
                                 *(c.get(r+'_'+s, 0) for r in ('guide', 'gene') for s in ('correct', 'incorrect', 'ambiguous', 'none')), c.get('invalid', 0)))
    write_tsv(output/'read-resolution.tsv', ['accession', 'policy', 'resolution', *STATES], qc_rows)
    write_tsv(output/'known-origin-controls.tsv', ['accession', 'policy', 'source_records', 'origin_in_candidates', 'guide_correct', 'guide_incorrect', 'guide_ambiguous', 'guide_none', 'gene_correct', 'gene_incorrect', 'gene_ambiguous', 'gene_none', 'invalid'], control_rows)
    detail_header = ['gene_annotation', 'policy', 'eligible_primary', 'fixed_exact_baseline_sum', 'baseline_guides_ge10', 'fixed_best_control', 'fixed_best_cellular', 'joint_control_lower', 'joint_cellular_lower', 'joint_control_upper', 'joint_cellular_upper', 'fixed_best_log2_ratio', 'joint_gene_unique_log2_ratio', 'delta_log2', 'conditional_low', 'conditional_high']
    write_tsv(output/'all-gene-representation.tsv', detail_header, detail)
    write_tsv(output/'all-primary-outliers.tsv', detail_header, sorted((r for r in detail if r[2] and abs(r[13]) >= 0.5), key=lambda r: (-abs(r[13]), r[0], r[1])))
    write_tsv(output/'effect-sensitivity.tsv', list(effects[0]), ([r[k] for k in effects[0]] for r in effects))
    summary = {'schema': 'dotmatch.research.ar003.report.v1', 'completion': 'complete',
               'audit_run': args.audit_run, 'audit_commit': AUDIT_COMMIT, 'protocol_commit': PROTOCOL_COMMIT,
               'full_public_records': sum(joint[a][1]['records'] for a in ACCESSIONS), 'archives': 3, 'assay_designs': 2,
               'biological_replicated_contrasts': 0, 'samples': summary_samples,
               'primary_effects': [r for r in effects if r['primary']], 'all_sensitivity_settings': effects,
               'validation': {'candidate_classes_reaggregated': class_counts, 'old_and_new_artifact_integrity': integrity,
                              'independent_SQL_Python_legacy_gene_aggregations': 2,
                              'ratio_vs_log_difference_arithmetic_checks': arithmetic,
                              'complete_class_to_guide_and_gene_reconciliation': True,
                              'source_sha256': digest(Path(__file__))},
               'limitations': ['Yusa plasmid/ESC1 is not a replicated biological contrast; Brunello is a plasmid sample.',
                              'Global multi-position and fixed-window comparisons also change extraction; this is not solely mismatch policy.',
                              'Additional gene-resolved records are not additional molecules or a demonstrated gain in accuracy.',
                              'Gene and guide views overlap and cannot be added. Gene labels are original annotations.',
                              'Bounds cover the candidate-supported assigned subset under the declared model, not all true source molecules.',
                              'Conditional ranges are not confidence intervals and their separate upper endpoints need not be jointly feasible.',
                              'Known-origin controls are constructed, error-free and based on previously inspected assay contexts.',
                              'No production release, manuscript submission, biological discovery or calibrated FDR claim.']}
    dump(output/'summary.json', summary)
    headline = []
    for accession in ACCESSIONS:
        s = summary_samples[accession]
        for p in POLICIES:
            headline.append([accession, p, number(s['qc'][p]['guide']['unique']), number(s['qc'][p]['gene']['unique']), number(s['additional_gene_resolved_records'][p])])
    primary_table = [[r['policy'], number(r['eligible_annotations']), number(r['shift_ge_half_log2']), f"{100*r['shift_fraction']:.3f}%", number(r['strong_sign_reversals'])] for r in summary['primary_effects']]
    report = '# Joint guide/position resolution in pooled CRISPR counting\n\n'
    report += '**AR003 technical results; 6 September 2026. Not peer reviewed.**\n\n'
    report += 'The new joint decoder was executed on all **30,215,791 original sequencing records** across ERR376998, ERR376999 and SRR8297997. It considered every permitted position jointly, rather than adding a count for each accepted window. Complete fixed-window baseline counts were independently reconciled with the prior pinned DotMatch 0.5.0 results.\n\n'
    report += '## What is identifiable at each resolution?\n\n'
    report += markdown_table(['Archive', 'Policy', 'Guide-unique reads', 'Gene-unique reads', 'Additional gene-identifiable reads'], headline)+'\n\n'
    report += 'The last column is gene-unique minus guide-unique within the same policy. These are overlapping views of the same records, not totals to add together. A gene-identifiable read may remain ambiguous between guides of the same gene. More identifiable reads do not, by themselves, establish more accurate assignments. Complete ambiguous, unmatched, invalid and position-level budgets are in `read-resolution.tsv`.\n\n'
    report += '## Yusa representation sensitivity\n\n'
    report += markdown_table(['Joint policy versus fixed best', 'Eligible annotation groups', 'Absolute change >=0.5 log2', 'Share', 'Strong sign reversals'], primary_table)+'\n\n'
    report += 'Eligibility is fixed-window exact plasmid sum >=50 with at least two guides each having >=10 baseline counts. The same population is used for every policy. Ratios use original input-read exposure and pseudocount 0.5. Strong sign reversals are a descriptive secondary check requiring opposite signs with both absolute ratios >=0.5. Thresholds 20/100 and pseudocounts 0.1/1 are reported separately. This changes extraction and resolution as well as assignment; it is not a pure mismatch-policy comparison.\n\n'
    report += f"All {len(new_genes['ERR376998']):,} original Yusa annotation labels are retained in the complete tables, including zero counts and ineligible groups. They are original annotation groups, not automatically distinct validated genes. `all-primary-outliers.tsv` contains every qualifying outlier, not only illustrative examples. No gene-level p-values or phenotype discoveries are inferred from this unreplicated pair.\n\n"
    report += '## Known-origin controls\n\n'
    control_table = [[r[0], r[1], number(r[2]), number(r[4]), number(r[5]), number(r[8]), number(r[9])] for r in control_rows]
    report += markdown_table(['Archive','Policy','Constructed records','Guide correct','Guide incorrect','Gene correct','Gene incorrect'], control_table)+'\n\n'
    report += 'These are the previously archived balanced, error-free constructs, one per reference guide. They test implementation and information loss under known contexts; they do not estimate real sequencing accuracy or generalize to unknown assays. Ambiguous and unassigned results are fully retained in `known-origin-controls.tsv`.\n\n'
    report += '## Validation and evidence boundaries\n\n'
    native_cells = sum(summary_samples[a]['validation']['full_fixed_window_native_count_cells'] for a in ACCESSIONS)
    oracle_records = sum(summary_samples[a]['validation']['exhaustive_full_library_full_position_records'] for a in ACCESSIONS)
    report += f'- **{native_cells:,}** complete fixed-window guide-count cells reconciled with pinned native evidence, plus all fixed-window read states.\n'
    report += f'- **{oracle_records:,}** fixed-seed selected public records checked against every reference target at every permitted position; these are sampled checks, not an exhaustive all-target check on every archive record.\n'
    report += f'- **{sum(class_counts.values()):,}** candidate classes independently reaggregated into every guide count, gene lower/upper bound and matched-state total.\n'
    report += f'- **{arithmetic:,}** representation calculations checked through ratio and log-difference arithmetic; old annotation sums independently reconciled through SQL and Python.\n\n'
    report += 'Full input SHA-256, MD5, byte count and record count were verified before execution and input SHA-256 checked again afterward. The report verifies each completed artifact checksum before analysis. Hash agreement proves byte identity, not biological accuracy. No unfinished local execution was used as evidence.\n\n'
    report += 'Conditional gene-count and effect ranges assume that each retained record originated from one of its candidates. They exclude unknown origins outside the model, including unmatched records. They are not confidence intervals, do not certify biological truth, and their marginal maxima need not be simultaneously achievable. No flanking-sequence validation, indel model, base-quality weighting, cell calling or UMI deduplication is implemented.\n\n'
    report += '## Scientific contribution and next publication gate\n\n'
    report += 'This establishes an executed read-conserving, joint-position audit with separate guide and gene resolution, not a new biological mechanism. The scientifically testable contribution is whether position-aware candidate accounting preserves useful gene-level information while exposing recycled guide evidence. Context-aware comparators, outcome-blind replicated independent screens and accuracy/calibration evidence remain necessary before claims of superior screening performance or a publication-ready general method.\n\n'
    report += 'Grouping ambiguous reads is established prior art; elementary marginal bounds are not claimed as new mathematics. The strongest potential paper combines an assay-specific mechanism, complete empirical accounting, validated implementation and downstream robustness evidence rather than merely observing that counting programs differ.\n\n'
    report += '## Reproduction\n\n'
    report += f'Protocol commit: `{PROTOCOL_COMMIT}`. Full-audit code commit: `{AUDIT_COMMIT}`. [Full audit run](https://github.com/dnncha/dotmatch/actions/runs/{args.audit_run}). [Prior independent replay](https://github.com/dnncha/dotmatch/actions/runs/34030217143). Source commands and hashes are recorded in every completion manifest. Primary source accessions remain the original input provenance; GitHub transport artifacts can expire.\n'
    (output/'REPORT.md').write_text(report, encoding='utf-8')
    validation = '# Validation assessment\n\n**Share with caveats: executed computational and descriptive results. Not ready for biological discovery or accuracy-superiority claims.**\n\n'
    validation += 'Every input artifact manifest, archive identity, guide identity, gene annotation join, candidate-class reconstruction, integer count budget and output bound was checked. The main report separates complete count reconciliation from sampled all-target validation. The unreplicated sample design and model-conditional nature of the ranges are material limitations, not footnotes.\n\n'
    validation += 'The local runtime failed during earlier full jobs. Those unfinished jobs were discarded as evidence. The reported full runs were independently executed in GitHub Actions from the recorded commit. The original standalone 127-test suite was rerun locally before failure; the integrated joint decoder has its own independently executed CI tests. These are different test scopes and are not summed into a fictitious single suite.\n\n'
    validation += 'Remaining blockers: independent replicated phenotype contrasts, comparable assay-aware decoders, stronger known-origin error models or orthogonal experimental validation, human review of manuscript and annotation interpretation. No claim that fixed-window counting is optimal for staggered Brunello reads.\n'
    (output/'VALIDATION.md').write_text(validation, encoding='utf-8')
    body = '<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>DotMatch AR003 research evidence</title><style>body{max-width:1100px;margin:48px auto;padding:0 24px;font:17px/1.6 system-ui,sans-serif}pre{white-space:pre-wrap;overflow-wrap:anywhere;font:15px/1.7 ui-monospace,monospace}h1{font-size:2rem}a{color:inherit}</style><h1>DotMatch AR003: complete research record</h1><p>Validated numerical report. Exploratory technical results, not peer reviewed.</p><pre>'+html.escape(report)+'</pre></html>\n'
    (output/'report.html').write_text(body, encoding='utf-8')
    dump(output/'completion.json', {'completion': 'complete', 'audit_run': args.audit_run,
                                    'files': {p.name: digest(p) for p in sorted(output.iterdir()) if p.is_file()}})
    print(json.dumps({'completion': 'complete', 'primary_effects': summary['primary_effects'],
                      'additional_gene_resolved_records': {a: summary_samples[a]['additional_gene_resolved_records'] for a in ACCESSIONS},
                      'validation': summary['validation']}, sort_keys=True), flush=True)


class Tests(unittest.TestCase):
    def test_log_formula_and_exposure(self):
        self.assertAlmostEqual(log_ratio(10, 10, 100, 100, .5), 0)
        self.assertAlmostEqual(log_ratio(10, 10, 200, 100, .5), 1)
        self.assertLess(log_ratio(0, 100, 1000, 1000, .5), 0)
        for args in ((-1, 2, 3, 4, .5), (1, 2, 0, 4, .5), (1, 2, 3, 4, 0)):
            with self.assertRaises(ValueError): log_ratio(*args)

    def test_strict_counts(self):
        self.assertEqual(natural('0'), 0)
        for value in ('-1', '1.5', 'nan', ' 2', '2 ', '\u0661', 1):
            with self.assertRaises(ValueError): natural(value)

    def test_duplicate_identifiers_rejected(self):
        with self.assertRaises(ValueError): keyed([{'id':'a'}, {'id':'a'}], 'id')

    def test_sql_aggregation(self):
        rows = [dict(id='a', gene='g', exact='10', radius_k1='9', best_k1='12'), dict(id='b', gene='g', exact='20', radius_k1='19', best_k1='22')]
        self.assertEqual(legacy_groups(rows)['g'], dict(exact=30, radius_k1=28, best_k1=34, supported_guides=2))

    def test_bounds_endpoints(self):
        nc, nt, pseudo = 100, 200, .5
        low = log_ratio(1, 20, nc, nt, pseudo)
        high = log_ratio(10, 2, nc, nt, pseudo)
        for c in range(2, 21):
            for t in range(1, 11):
                point = log_ratio(t, c, nc, nt, pseudo)
                self.assertLessEqual(low, point)
                self.assertLessEqual(point, high)

    def test_markdown_control_table(self):
        rendered = markdown_table(['Archive','Correct','Incorrect'], [['ERR', 10, 0], ['SRR', 20, 0]])
        self.assertIn('| ERR | 10 | 0 |', rendered)
        self.assertEqual(len(rendered.splitlines()), 4)


def main():
    if '--test' in sys.argv:
        unittest.main(argv=[sys.argv[0]], verbosity=2)
        return
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--joint', required=True)
    p.add_argument('--replay', required=True)
    p.add_argument('--out', required=True)
    p.add_argument('--audit-run', type=int, default=AUDIT_RUN)
    build(p.parse_args())


if __name__ == '__main__':
    main()
