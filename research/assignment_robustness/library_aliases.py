#!/usr/bin/env python3
"""Post-hoc exact-overlap library audit, declared before these graph results.

Motivated by observed Yusa multi-offset records mostly sharing gene annotation.
For shifts 1, 2 and 3, enumerate ordered guide pairs A,B with A[s:] == B[:-s].
This is a structural descriptor, NOT a complete Hamming-one ambiguity graph,
a predictor validated on held-out data, or proof of biological misassignment.
Flank sequence compatibility and errors can add or remove observed aliases.
"""
from __future__ import annotations
import argparse
from collections import defaultdict
import itertools
import json
from pathlib import Path
import tempfile
import unittest
import pilot


def aliases(library, shifts=(1, 2, 3)):
    rows = []
    length = len(library[0][1])
    pilot.require(all(len(row[1]) == length for row in library), 'mixed reference lengths')
    for shift in shifts:
        pilot.require(0 < shift < length, 'invalid overlap shift')
        prefixes = defaultdict(list)
        for i, (_ident, seq, _gene) in enumerate(library):
            prefixes[seq[:-shift]].append(i)
        for a, (ident_a, seq, gene_a) in enumerate(library):
            for b in prefixes.get(seq[shift:], []):
                ident_b, _seq_b, gene_b = library[b]
                rows.append((shift, ident_a, ident_b, gene_a, gene_b, int(gene_a == gene_b), int(a == b)))
    return rows


def run(out):
    out.mkdir(parents=True, exist_ok=False)
    work = out / 'work'
    work.mkdir()
    summary = {'schema': 'dotmatch.research.exact_overlap.v1', 'status': 'running',
               'analysis_kind': 'post-hoc structural explanation, not held-out predictive validation',
               'harness_sha256': pilot.digest(__file__), 'shifts': [1, 2, 3], 'libraries': {}}
    try:
        for name in ('yusa', 'brunello'):
            library, path, source = pilot.fetch_library(name, work)
            rows = aliases(library)
            pilot.table(out / f'{name}.exact_overlap_edges.tsv',
                ['shift', 'guide_a_at_earlier_offset', 'guide_b_at_later_offset', 'gene_a', 'gene_b', 'same_gene', 'same_guide'], rows)
            summary['libraries'][name] = {'source': source, 'guide_rows': len(library),
                'ordered_edges': len(rows), 'cross_gene_edges': sum(not row[5] for row in rows),
                'same_guide_edges': sum(row[6] for row in rows),
                'guides_in_edges': len({guide for row in rows for guide in row[1:3]}),
                'edges_by_shift': {str(s): sum(row[0] == s for row in rows) for s in (1, 2, 3)}}
        summary['status'] = 'complete'
    except Exception as exc:
        summary['status'] = 'failed'
        summary['error'] = f'{type(exc).__name__}: {exc}'
        raise
    finally:
        pilot.save_json(out / 'run.json', summary)
        pilot.save_json(out / 'MANIFEST.json', {'status': summary['status'], 'files': {
            p.name: {'sha256': pilot.digest(p), 'bytes': p.stat().st_size}
            for p in out.iterdir() if p.is_file() and p.name != 'MANIFEST.json'}})
        print(json.dumps(summary), flush=True)


class Tests(unittest.TestCase):
    def test_exhaustive_graph(self):
        library = [(str(i), ''.join(s), str(i % 5)) for i, s in enumerate(itertools.product('ACGT', repeat=4))]
        observed = {(row[0], row[1], row[2]) for row in aliases(library)}
        expected = {(shift, a[0], b[0]) for shift in (1, 2, 3) for a in library for b in library if a[1][shift:] == b[1][:-shift]}
        self.assertEqual(observed, expected)

    def test_self_and_duplicate_sequences_preserved(self):
        library = [('a', 'AAAA', 'G'), ('b', 'AAAA', 'H')]
        self.assertEqual(len(aliases(library)), 12)
        self.assertEqual(sum(row[-1] for row in aliases(library)), 6)


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == 'test':
        unittest.main(argv=[sys.argv[0]], verbosity=2)
    else:
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument('--out', required=True)
        run(Path(parser.parse_args().out).resolve())
