#!/usr/bin/env python3
"""Validate three complete audit packets and publish only aggregate evidence."""
from __future__ import annotations
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()


def save(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')


def main(source, destination):
    source, destination = Path(source), Path(destination)
    records = {}
    for completion in source.rglob('completion.json'):
        if completion.parent.name != 'audit-output':
            continue
        row = json.loads(completion.read_text())
        if row['completion'] != 'complete' or row['assignment_disagreements'] != 0:
            raise ValueError('Incomplete or non-reconciled input')
        if row['accession'] in records:
            raise ValueError('Duplicate completed sample')
        for name, expected in row['files'].items():
            p = completion.parent / name
            if not p.resolve().is_relative_to(completion.parent.resolve()) or digest(p) != expected:
                raise ValueError('Artifact member path/digest mismatch')
        if row['statistics']['reads'] != int(row['input']['read_count']):
            raise ValueError('Record conservation differs')
        records[row['accession']] = (row, completion.parent)
    order = ['ERR376998', 'ERR376999', 'SRR8297997']
    if set(records) != set(order):
        raise ValueError('Missing one of three complete archives')
    expected_yusa = {
        'ERR376998': (10093905, 10161871, 756973, 720995),
        'ERR376999': (10300758, 10108301, 742411, 707092),
    }
    for accession, expected in expected_yusa.items():
        s = records[accession][0]['statistics']
        observed = tuple(s.get(k, 0) for k in ('reads', 'count_events', 'extra_events', 'multiple_distinct_guides'))
        if observed != expected:
            raise ValueError(f'Previous Yusa result did not reproduce: {accession}: {observed}')
    destination.mkdir(parents=True, exist_ok=True)
    selected = [records[a][0] for a in order]
    run = os.environ['GITHUB_RUN_ID']
    run_url = f'https://github.com/dnncha/dotmatch/actions/runs/{run}'
    totals = {k: sum(r['statistics'].get(k, 0) for r in selected) for k in
              ('reads', 'matched_reads', 'count_events', 'extra_events', 'multiple_distinct_guides',
               'within_gene_multiple', 'cross_gene_multiple', 'repeated_same_target_events')}
    yusa = {k: sum(r['statistics'].get(k, 0) for r in selected[:2]) for k in totals}
    for row, folder in records.values():
        short = dict(row)
        short.pop('files')
        # Full per-file provenance and tables remain in the downloadable job artifact.
        save(destination / (row['accession'] + '.json'), short)
        shutil.copy2(folder / 'execution-start.json', destination / (row['accession'] + '.execution.json'))
    aggregate = {'completion': 'complete', 'scope': 'three_complete_archives_two_library_constructs',
                 'totals': totals, 'yusa_totals': yusa, 'workflow_run': run,
                 'workflow_url': run_url, 'biological_replicate_count_claim': False,
                 'all_target_oracle_windows': sum(r['all_target_oracle_windows'] for r in selected),
                 'all_guide_native_comparisons': sum(r['all_guide_native_comparisons'] for r in selected),
                 'all_guide_upstream_comparisons': sum(r['all_guide_upstream_comparisons'] for r in selected),
                 'assignment_disagreements': 0,
                 'source_sha256': {name: digest(Path(__file__).parent / name)
                                   for name in ('replay.py', 'publish_results.py', 'REPLAY_PROTOCOL.md')}}
    save(destination / 'aggregate.json', aggregate)
    columns = ['accession', 'reference_targets', 'reads', 'matched_reads', 'count_events', 'extra_events',
               'multiple_distinct_guides', 'within_gene_multiple', 'cross_gene_multiple']
    with (destination / 'read-conservation.tsv').open('w', newline='') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n'); w.writerow(columns)
        for r in selected:
            w.writerow([r['accession'], r['target_count'], *(r['statistics'].get(k, 0) for k in columns[2:])])
    save(destination / 'report-rendering.json', {
        'selected_audience': 'technical', 'selected_report_mode': 'html',
        'status': 'blocked', 'reason': 'The standard analytics portable-report builder was not accessible after the local execution service failed. No custom substitute is misrepresented as that renderer.',
        'delivered_alternative': 'Audited Markdown research report, source code, aggregate TSV/JSON, and complete workflow artifacts.'})
    lines = ['# One read, several guides', '', '## Technical summary', '',
        f'**Executed evidence, 6 September 2026:** independent full-read audits of {totals["reads"]:,} reads across three public archives and two library/construct designs. Every fixed-window guide count reconciled with pinned DotMatch 0.5.0; every multi-offset count reconciled with the unchanged guide-counter v0.1.3 executable. No assignment disagreements were found.', '',
        f'In the two Yusa archives, {yusa["multiple_distinct_guides"]:,} reads entered more than one guide row. Of these, {yusa["within_gene_multiple"]:,} ({100*yusa["within_gene_multiple"]/yusa["multiple_distinct_guides"]:.2f}%) involved guides labelled with the same gene. There were {yusa["extra_events"]:,} extra count events. These are computational read-reuse measurements, not estimates of incorrect biological origins.', '',
        '**Interpretation:** apparent support from several guide rows can reuse the same reads. Fixed-position Hamming separation alone does not diagnose aliases created by shifted extraction windows. This is an assay-dependent effect, not evidence that every screen or every guide-counting program is affected.', '',
        '## Scope and definitions', '',
        'ERR376998 and ERR376999 are the Yusa plasmid/ESC1 tutorial pair, not independent replicated treatment arms. SRR8297997 is a Brunello plasmid library run, not a treatment contrast. The analysis counts sequencing records, not deduplicated original molecules; no UMI deduplication or gene-level hypothesis testing is performed.', '',
        'Let N be input reads, M reads with at least one accepted window, C retained guide-window count events, U=N-M unmatched reads, and E=C-M extra events. Then C-N=E-U. A count-table total below the input-read total does not rule out repeated counting. A multiple-guide read is a record with at least two distinct accepted target IDs; it is not necessarily a cross-gene ambiguity.', '',
        '## Complete-archive results', '',
        '| Archive | Reads N | Matched reads M | Count events C | Extra events E | Multiple-guide reads | Same-gene cases | Cross-gene cases |',
        '|---|---:|---:|---:|---:|---:|---:|---:|']
    for r in selected:
        s = r['statistics']
        lines.append('| ' + r['accession'] + ' | ' + ' | '.join(f'{s.get(k,0):,}' for k in
                     ('reads','matched_reads','count_events','extra_events','multiple_distinct_guides','within_gene_multiple','cross_gene_multiple')) + ' |')
    lines += ['', f'For Yusa, the pooled count-event total is {yusa["count_events"]:,}, below {yusa["reads"]:,} input reads, despite {yusa["extra_events"]:,} excess events. Unmatched reads conceal the duplication in this aggregate check.', '',
        'The table is a census of the retrieved archives under the specified rules. The sample-to-sample comparison is descriptive; treating millions of reads as millions of biological replicates would be invalid.', '',
        '## Discovery-only construct model and held-out evaluation', '',
        'Each sample uses its first 100,000 records solely for offset discovery and modal-flank estimation. Tied flank bases become N. Each reference guide is inserted into that inferred context, and its accepted guide-window events are recorded. The entire model is written before later records are evaluated. The prediction target is whether a read produces multiple count events—not its true biological source.', '',
        '| Archive | Templates predicting multiple events | Evaluated exact-window reads | Excluded evaluation reads | TP | FP | FN | TN |',
        '|---|---:|---:|---:|---:|---:|---:|---:|']
    for r in selected:
        p = r['template_prediction']; c = p['confusion']
        lines.append('| ' + r['accession'] + ' | ' + ' | '.join(f'{v:,}' for v in
                     [p['predicted_multiple_templates'],p['evaluation_reads'],p['excluded_evaluation_reads'],*(c.get(k,0) for k in ('TP','FP','FN','TN'))]) + ' |')
    lines += ['', 'The model is evaluated only where the chosen fixed window has a unique exact target. Excluded records remain in the full-read audit. This restriction and the near-constant construct sequence make the mechanistic prediction task much easier than biological-origin inference; no general assignment-accuracy claim follows.', '',
        '## Error-free interventions isolate the mechanism', '',
        'The following controls contain constructed reads with known source guide labels. They are not newly measured biological samples. Balanced controls contain exactly one error-free template per reference guide. The Yusa PHF23 source-only witness was selected after inspecting the real-data findings and is explicitly a post-hoc mechanistic test.', '',
        '| Reference/sample | Synthetic experiment | Upstream rule | Known source reads | Count events | Zero-origin guides receiving counts |',
        '|---|---|---|---:|---:|---:|']
    for r in selected:
        for x in r['synthetic_controls']:
            lines.append(f'| {r["accession"]} | {x["experiment"]} | {x["policy"]} | {x["source_reads"]:,} | {x["count_events"]:,} | {x["zero_origin_guides_receiving_counts"]:,} |')
    lines += ['', 'In the selected witness, the two 19-base targets differ at 17 positions when compared in the same frame. Nevertheless the construct flank plus an 18-base target overlap yields an exact match to the second guide at a shifted offset. Changing the last left-flank base removes that exact alias while a one-mismatch match can remain. Fixed-window exact DotMatch counts agree with the known source counts in every synthetic control.', '',
        'An eightfold change in source-only synthetic reads produces the same apparent change in the alias guide although that guide contributes no source reads. This is a demonstration of duplicated guide-level evidence, not a new gene-function finding or a calibrated false-discovery-rate result.', '',
        '## Methods and independent checks', '',
        'The full input archives were verified against ENA byte counts and MD5s, then locked with SHA-256. The corrected references retain target and gene identities. Fixed-window policy comparisons use zero-based start/length 23/19 for Yusa and 21/20 for Brunello. These are controlled policy-comparison windows, not claims of optimal extraction for staggered reads.', '',
        'The independent index uses two disjoint exact seeds: any observation within one Hamming substitution must share at least one complete seed. Every candidate is checked by full Hamming distance, with duplicate target identities retained. It does not use DotMatch candidate lists. A 12,500-window constructed grid checks all candidates and all three policies against exhaustive all-target enumeration, including N and duplicate cases. The real-read audit also checks 200 prespecified pseudorandom ordinals per archive against every reference target.', '',
        f'The replay reconciles {aggregate["all_guide_native_comparisons"]:,} native policy/guide count values and {aggregate["all_guide_upstream_comparisons"]:,} upstream guide count values, with {aggregate["all_target_oracle_windows"]:,} exhaustive real-window checks. Each value is derived from the complete archive, not just the oracle sample. Count-table identity, gene labels, all four read states, changed-read totals, and event-conservation identities are checked.', '',
        'Upstream reproduction intentionally includes its exact-match precedence, last-inserted duplicate-sequence behaviour, ACGT-only mismatch lookup, event-based offset threshold and increment-at-every-accepted-offset loop. These semantics are not silently equated with fixed-window Hamming matching or with true molecular origin.', '',
        '## Prior art and what could be publishable', '',
        'Mapping uncertainty is not a new idea: [bcSeq](https://doi.org/10.1093/bioinformatics/bty402) models sequencing-error-aware barcode assignment, and [crispat](https://doi.org/10.1093/bioinformatics/btae535) studies downstream sensitivity to guide-to-cell assignment. [Buschmann and Bystrykh](https://doi.org/10.1186/1471-2105-14-272) already describe why DNA sequence context matters for barcode decoding. The general boundary/uncertainty concepts must not be claimed as novel.', '',
        '**The strongest candidate contribution is a bulk-screen, read-to-count evidence audit showing when apparent multi-guide support reuses the same records, an assay-context diagnostic, and independent validation that separates this effect from sequencing errors and mismatch rescue.** The current study supports a methods/technical-results draft. A stronger downstream-statistics paper still requires replicated biological contrasts, a locked normalization/filtering plan, held-out studies, and direct tests of hit/rank/effect changes.', '',
        '## Limitations and unresolved questions', '',
        'Two library/construct designs are not an atlas of all CRISPR screens. The Yusa results motivated part of the subsequent mechanism work; that adaptation is documented rather than presented as preregistered discovery. A single-read witness cannot establish empirical prevalence, and algorithm agreement cannot certify biological truth. Guide-window aliases can be same-gene without being harmless to guide-level support, but gene-level statistical distortion has not been established here.', '',
        'The independently generated original pilot bounds remain conditional on candidate-set completeness. They are not confidence intervals and are not used to assert calibrated gene-level coverage in this replay. Likewise, no confidence intervals from millions of read records are substituted for biological replication.', '',
        'The next decisive question is whether the documented reuse changes replicated guide-support and gene-level conclusions under otherwise identical analysis. Full differential-screen curation, outcome-blind holdouts and appropriate comparator coverage remain publication gates, not completed claims.', '',
        '## Reproducibility and evidence access', '',
        f'[Completed workflow and complete per-sample artifacts]({run_url}). Source and protocols live beside this report. Job artifacts contain complete count tables, template models, read-class tables, synthetic controls, commands, source and binary hashes, and failure-visible logs. Public raw FASTQs are not recommitted. Archive transport artifacts expire; original ENA accessions and acquisition scripts remain the source of input data.', '',
        'The formal portable HTML renderer was unavailable after a local execution-service failure. This Markdown report is the explicitly labelled alternative, not a claim of completed HTML rendering. No production release, manuscript submission, accepted paper or new biological discovery is claimed.']
    (destination / 'README.md').write_text('\n'.join(lines) + '\n')
    save(destination / 'checksums.json', {p.name: digest(p) for p in sorted(destination.iterdir()) if p.is_file() and p.name != 'checksums.json'})
    print(json.dumps(aggregate, indent=2))

if __name__ == '__main__':
    if len(sys.argv) != 3:
        raise SystemExit('usage: publish_results.py COMPLETED_ARTIFACTS OUTPUT_DIRECTORY')
    main(*sys.argv[1:])
