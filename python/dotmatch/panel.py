from __future__ import annotations

import argparse
import csv
import html
import json
import math
import random
import time
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from itertools import combinations
from pathlib import Path
from typing import Iterable, Iterator, Sequence


DNA = "ACGT"
COMPLEMENT = str.maketrans("ACGT", "TGCA")
PANEL_SCHEMA_VERSION = 1
MAX_EXACT_CERTIFICATE_RADIUS = 2


@dataclass(frozen=True)
class PanelRecord:
    barcode_id: str
    sequence: str


@dataclass(frozen=True)
class AssignmentOutcome:
    target_index: int
    best_distance: int
    second_best_distance: int
    match_count: int
    status: str


@dataclass
class AssignmentConfig:
    metric: str = "hamming"
    k: int = 1
    ambiguity_policy: str = "radius"
    fixed_window: bool = True
    reverse_complement_mode: str = "warn"


@dataclass
class PanelConstraints:
    length: int | None = None
    alphabet: str = DNA
    allow_ambiguous_literals: bool = False
    min_hamming_distance: int = 5
    min_levenshtein_distance: int = 4
    min_sequence_levenshtein_distance: int = 4
    gc_min: float = 0.35
    gc_max: float = 0.65
    max_homopolymer: int = 3
    max_dinucleotide_repeat: int = 4
    avoid_reverse_complements: bool = True
    avoid_self_complementarity: bool = True
    avoid_prefix_collisions: bool = True
    avoid_suffix_collisions: bool = True
    avoid_ambiguous_bases: bool = True
    forbidden_motifs: list[str] = field(default_factory=lambda: ["AAAA", "CCCC", "GGGG", "TTTT", "GATC", "AAGCTT"])
    left_flank: str = ""
    right_flank: str = ""
    context_window: int = 0
    cycle_balance_enabled: bool = True
    cycle_min_base_fraction: float = 0.15
    cycle_max_base_fraction: float = 0.40


@dataclass
class DesignOptions:
    panel_name: str = "dotmatch_panel"
    count: int = 96
    length: int = 16
    seed: int = 1
    engine: str = "greedy"
    restarts: int = 100
    candidate_pool_size: int = 100000
    objective: str = "strict-demux"
    max_runtime: float | None = None
    target_count: int | None = None
    population: int = 500
    iterations: int = 1000
    mutation_rate: float = 0.05
    elite_fraction: float = 0.10
    experimental_signal: bool = False
    chemistry: str = ""
    plate_format: int = 96


PRESETS: dict[str, dict[str, object]] = {
    "strict-24x12": {
        "count": 24,
        "length": 12,
        "min_hamming_distance": 5,
        "min_levenshtein_distance": 4,
        "min_sequence_levenshtein_distance": 4,
        "objective": "strict-demux",
        "description": "small panels; very high spacing",
    },
    "strict-96x16": {
        "count": 96,
        "length": 16,
        "min_hamming_distance": 5,
        "min_levenshtein_distance": 4,
        "min_sequence_levenshtein_distance": 4,
        "objective": "strict-demux",
        "description": "general inline barcodes; safe one-edit correction",
    },
    "strict-384x20": {
        "count": 384,
        "length": 20,
        "min_hamming_distance": 5,
        "min_levenshtein_distance": 4,
        "min_sequence_levenshtein_distance": 4,
        "objective": "max-count",
        "description": "large panels; high capacity",
    },
    "illumina-inline-96": {
        "count": 96,
        "length": 10,
        "min_hamming_distance": 4,
        "min_levenshtein_distance": 3,
        "min_sequence_levenshtein_distance": 3,
        "objective": "strict-demux",
        "description": "fixed-window inline demux; Hamming-first",
    },
    "illumina-inline-strict": {
        "count": 96,
        "length": 16,
        "min_hamming_distance": 5,
        "min_levenshtein_distance": 4,
        "min_sequence_levenshtein_distance": 4,
        "objective": "strict-demux",
        "description": "strict fixed-window inline demux",
    },
    "illumina-dual-384": {
        "count": 384,
        "length": 10,
        "min_hamming_distance": 4,
        "min_levenshtein_distance": 4,
        "min_sequence_levenshtein_distance": 4,
        "objective": "strict-demux",
        "description": "paired i7/i5 sample indexing; pair-aware safety",
    },
    "nanopore-indel-robust-24": {
        "count": 24,
        "length": 18,
        "min_hamming_distance": 5,
        "min_levenshtein_distance": 5,
        "min_sequence_levenshtein_distance": 5,
        "objective": "nanopore-indel-robust",
        "description": "long-read noisy contexts; Levenshtein/seqlev spacing",
    },
    "ont-rna004-signal-12": {
        "count": 12,
        "length": 16,
        "min_hamming_distance": 5,
        "min_levenshtein_distance": 5,
        "min_sequence_levenshtein_distance": 5,
        "objective": "ont-signal-separable",
        "description": "experimental symbolic + signal margin candidate narrowing",
    },
    "feature-barcode-192": {
        "count": 192,
        "length": 16,
        "min_hamming_distance": 5,
        "min_levenshtein_distance": 4,
        "min_sequence_levenshtein_distance": 4,
        "objective": "strict-demux",
        "description": "feature barcode fixed-window panels",
    },
    "crispr-guide-tag-96": {
        "count": 96,
        "length": 12,
        "min_hamming_distance": 5,
        "min_levenshtein_distance": 4,
        "min_sequence_levenshtein_distance": 4,
        "objective": "strict-demux",
        "description": "CRISPR guide tag panels",
    },
}


OBJECTIVE_WEIGHTS: dict[str, dict[str, float]] = {
    "strict-demux": {
        "min_pairwise_hamming_distance": 10.0,
        "min_pairwise_levenshtein_distance": 8.0,
        "min_pairwise_sequence_levenshtein_distance": 6.0,
        "pairwise_distance_5th_percentile": 4.0,
        "edit_sphere_separation_score": 4.0,
        "cycle_balance_score": 3.0,
        "plate_layout_score": 3.0,
        "gc_distribution_score": 2.0,
        "context_separation_score": 2.0,
        "ambiguous_variant_count": -8.0,
        "reverse_complement_collision_count": -6.0,
        "prefix_collision_count": -5.0,
        "forbidden_motif_count": -4.0,
        "low_complexity_penalty": -3.0,
    },
    "max-count": {
        "min_pairwise_hamming_distance": 7.0,
        "min_pairwise_levenshtein_distance": 5.0,
        "min_pairwise_sequence_levenshtein_distance": 4.0,
        "pairwise_distance_5th_percentile": 3.0,
        "edit_sphere_separation_score": 2.0,
        "cycle_balance_score": 2.0,
        "plate_layout_score": 1.0,
        "gc_distribution_score": 2.0,
        "context_separation_score": 1.0,
        "ambiguous_variant_count": -8.0,
        "reverse_complement_collision_count": -4.0,
        "prefix_collision_count": -4.0,
        "forbidden_motif_count": -4.0,
        "low_complexity_penalty": -3.0,
    },
    "nanopore-indel-robust": {
        "min_pairwise_hamming_distance": 6.0,
        "min_pairwise_levenshtein_distance": 11.0,
        "min_pairwise_sequence_levenshtein_distance": 10.0,
        "pairwise_distance_5th_percentile": 5.0,
        "edit_sphere_separation_score": 5.0,
        "cycle_balance_score": 2.0,
        "plate_layout_score": 1.0,
        "gc_distribution_score": 2.0,
        "context_separation_score": 4.0,
        "ambiguous_variant_count": -10.0,
        "reverse_complement_collision_count": -4.0,
        "prefix_collision_count": -6.0,
        "forbidden_motif_count": -4.0,
        "low_complexity_penalty": -3.0,
    },
    "ont-signal-separable": {
        "min_pairwise_hamming_distance": 6.0,
        "min_pairwise_levenshtein_distance": 8.0,
        "min_pairwise_sequence_levenshtein_distance": 8.0,
        "pairwise_distance_5th_percentile": 3.0,
        "edit_sphere_separation_score": 4.0,
        "cycle_balance_score": 2.0,
        "plate_layout_score": 1.0,
        "gc_distribution_score": 2.0,
        "context_separation_score": 4.0,
        "ambiguous_variant_count": -10.0,
        "reverse_complement_collision_count": -6.0,
        "prefix_collision_count": -5.0,
        "forbidden_motif_count": -4.0,
        "low_complexity_penalty": -3.0,
    },
}


def command_panel_namespace(argv: Sequence[str]) -> int:
    parser = argparse.ArgumentParser(
        prog="dotmatch panel",
        description="Design, certify, simulate, lay out, and export barcode panels for DotMatch known-target assignment.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    check = sub.add_parser("check", help="certify a barcode panel under configured DotMatch assignment semantics")
    check.add_argument("panel")
    _add_assignment_args(check)
    _add_constraint_args(check)
    _add_context_args(check)
    check.add_argument("--out-dir", required=True)

    design = sub.add_parser("design", help="design a new single-index barcode panel")
    design.add_argument("--spec")
    design.add_argument("--preset", choices=sorted(PRESETS))
    design.add_argument("--n", "--count", dest="count", type=int)
    design.add_argument("--length", type=int)
    design.add_argument("--seed", type=int)
    design.add_argument("--engine", choices=["greedy", "graph", "evolve"], default=None)
    design.add_argument("--restarts", type=int)
    design.add_argument("--candidate-pool-size", type=int)
    design.add_argument("--max-runtime", type=float)
    design.add_argument("--target-count", type=int)
    design.add_argument("--population", type=int)
    design.add_argument("--iterations", type=int)
    design.add_argument("--mutation-rate", type=float)
    design.add_argument("--elite-fraction", type=float)
    design.add_argument("--objective", choices=sorted(OBJECTIVE_WEIGHTS))
    design.add_argument("--experimental", action="store_true")
    design.add_argument("--experimental-signal", action="store_true")
    design.add_argument("--chemistry", default="")
    design.add_argument("--out-dir", required=True)
    _add_assignment_args(design)
    _add_constraint_args(design)
    _add_context_args(design)

    optimize = sub.add_parser("optimize", help="choose the safest subset from candidate/vendor barcodes")
    optimize.add_argument("candidates")
    optimize.add_argument("--n", "--count", dest="count", type=int, required=True)
    optimize.add_argument("--seed", type=int, default=1)
    optimize.add_argument("--engine", choices=["greedy", "graph", "evolve"], default="greedy")
    optimize.add_argument("--restarts", type=int, default=20)
    optimize.add_argument("--objective", choices=sorted(OBJECTIVE_WEIGHTS), default="strict-demux")
    optimize.add_argument("--out-dir", required=True)
    _add_assignment_args(optimize)
    _add_constraint_args(optimize)
    _add_context_args(optimize)

    simulate = sub.add_parser("simulate", help="stress-test assignment behavior under simple error models")
    simulate.add_argument("panel")
    simulate.add_argument("--reads", type=int, default=100000)
    simulate.add_argument("--seed", type=int, default=1)
    simulate.add_argument("--substitution-rate", type=float, default=0.005)
    simulate.add_argument("--insertion-rate", type=float, default=0.001)
    simulate.add_argument("--deletion-rate", type=float, default=0.001)
    simulate.add_argument("--quality-model", default="simple")
    simulate.add_argument("--out-dir", required=True)
    _add_assignment_args(simulate)

    layout = sub.add_parser("layout", help="create a plate-aware barcode layout")
    layout.add_argument("panel")
    layout.add_argument("--plate", type=int, choices=[24, 96, 384], default=96)
    layout.add_argument("--maximize-neighbor-distance", action="store_true", default=True)
    layout.add_argument("--avoid-row-column-near-neighbors", action="store_true", default=True)
    layout.add_argument("--out", required=True)
    _add_assignment_args(layout)

    export = sub.add_parser("export", help="write lab-ready sample-sheet files")
    export.add_argument("panel")
    export.add_argument("--format", choices=["illumina-samplesheet", "generic-tsv"], default="illumina-samplesheet")
    export.add_argument("--out-dir", required=True)
    _add_assignment_args(export)

    compare = sub.add_parser("compare", help="compare two barcode panels")
    compare.add_argument("old_panel")
    compare.add_argument("new_panel")
    compare.add_argument("--out-dir", required=True)
    _add_assignment_args(compare)

    dual = sub.add_parser("design-dual", help="design a paired i7/i5 barcode panel")
    dual.add_argument("--samples", type=int, required=True)
    dual.add_argument("--i7-count", type=int, required=True)
    dual.add_argument("--i5-count", type=int, required=True)
    dual.add_argument("--i7-length", type=int, required=True)
    dual.add_argument("--i5-length", type=int, required=True)
    dual.add_argument("--unique-dual", action="store_true")
    dual.add_argument("--combinatorial", action="store_true")
    dual.add_argument("--min-i7-distance", type=int, default=4)
    dual.add_argument("--min-i5-distance", type=int, default=4)
    dual.add_argument("--min-pair-distance", type=int, default=8)
    dual.add_argument("--seed", type=int, default=1)
    dual.add_argument("--out-dir", required=True)
    _add_assignment_args(dual)

    args = parser.parse_args(list(argv))
    try:
        if args.command == "check":
            return command_panel_check(args)
        if args.command == "design":
            return command_panel_design(args)
        if args.command == "optimize":
            return command_panel_optimize(args)
        if args.command == "simulate":
            return command_panel_simulate(args)
        if args.command == "layout":
            return command_panel_layout(args)
        if args.command == "export":
            return command_panel_export(args)
        if args.command == "compare":
            return command_panel_compare(args)
        if args.command == "design-dual":
            return command_panel_design_dual(args)
    except BrokenPipeError:
        return 1
    except Exception as exc:
        print(f"dotmatch panel: {exc}", file=__import__("sys").stderr)
        return 2
    parser.error("unreachable")
    return 2


def _add_assignment_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--k", type=int, default=1)
    parser.add_argument("--metric", choices=["hamming", "levenshtein", "sequence-levenshtein", "seqlev"], default="hamming")
    parser.add_argument("--reverse-complement-mode", choices=["ignore", "warn", "fail"], default="warn")


def _add_constraint_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--min-hamming-distance", type=int)
    parser.add_argument("--min-levenshtein-distance", type=int)
    parser.add_argument("--min-sequence-levenshtein-distance", type=int)
    parser.add_argument("--gc-min", type=float)
    parser.add_argument("--gc-max", type=float)
    parser.add_argument("--max-homopolymer", type=int)
    parser.add_argument("--max-dinucleotide-repeat", type=int)
    parser.add_argument("--avoid-rc", "--avoid-reverse-complements", dest="avoid_reverse_complements", action="store_true", default=None)
    parser.add_argument("--allow-ambiguous-literals", action="store_true")
    parser.add_argument("--forbidden-motif", action="append", dest="forbidden_motifs")


def _add_context_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--left-flank", default="")
    parser.add_argument("--right-flank", default="")
    parser.add_argument("--context-window", type=int, default=None)


def command_panel_check(args: argparse.Namespace) -> int:
    records = read_panel(args.panel)
    assignment = _assignment_from_args(args)
    constraints = _constraints_from_args(args, length=_single_length(records))
    out_dir = Path(args.out_dir)
    certificate = check_panel(records, assignment, constraints, out_dir, source_panel=Path(args.panel))
    return 0 if certificate["status"] != "fail" else 1


def command_panel_design(args: argparse.Namespace) -> int:
    options, assignment, constraints = _design_from_args(args)
    if options.objective == "ont-signal-separable" and not (options.experimental_signal or getattr(args, "experimental", False)):
        raise ValueError("--objective ont-signal-separable requires --experimental or --experimental-signal")
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    records, trace = design_panel(options, assignment, constraints)
    if len(records) < options.count:
        raise ValueError(f"could only design {len(records)} compatible barcodes; requested {options.count}")
    _write_design_outputs(out_dir, records, trace, options, assignment, constraints, mode="design")
    return 0


def command_panel_optimize(args: argparse.Namespace) -> int:
    candidates = read_panel(args.candidates)
    if not candidates:
        raise ValueError("candidate panel is empty")
    length = _single_length(candidates)
    assignment = _assignment_from_args(args)
    constraints = _constraints_from_args(args, length=length)
    options = DesignOptions(
        panel_name="optimized_panel",
        count=args.count,
        length=length or len(candidates[0].sequence),
        seed=args.seed,
        engine=args.engine,
        restarts=args.restarts,
        candidate_pool_size=len(candidates),
        objective=args.objective,
    )
    ranked = [record for record in candidates if candidate_filter_reason(record.sequence, constraints) == "pass"]
    if len(ranked) < args.count:
        ranked = candidates
    records, trace = select_panel_from_candidates(ranked, options, constraints)
    if len(records) < args.count:
        raise ValueError(f"could only optimize {len(records)} compatible barcodes; requested {args.count}")
    _write_design_outputs(Path(args.out_dir), records, trace, options, assignment, constraints, mode="optimize")
    return 0


def command_panel_simulate(args: argparse.Namespace) -> int:
    records = read_panel(args.panel)
    assignment = _assignment_from_args(args)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    simulate_panel(
        records,
        assignment,
        reads=args.reads,
        seed=args.seed,
        substitution_rate=args.substitution_rate,
        insertion_rate=args.insertion_rate,
        deletion_rate=args.deletion_rate,
        quality_model=args.quality_model,
        out_dir=out_dir,
        source_panel=Path(args.panel),
    )
    return 0


def command_panel_layout(args: argparse.Namespace) -> int:
    records = read_panel(args.panel)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    assignment = _assignment_from_args(args)
    write_plate_layout(records, assignment, args.plate, out)
    return 0


def command_panel_export(args: argparse.Namespace) -> int:
    records = read_panel(args.panel)
    assignment = _assignment_from_args(args)
    export_panel(records, assignment, Path(args.out_dir), fmt=args.format, source_panel=Path(args.panel))
    return 0


def command_panel_compare(args: argparse.Namespace) -> int:
    old_records = read_panel(args.old_panel)
    new_records = read_panel(args.new_panel)
    assignment = _assignment_from_args(args)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    compare_panels(old_records, new_records, assignment, out_dir, Path(args.old_panel), Path(args.new_panel))
    return 0


def command_panel_design_dual(args: argparse.Namespace) -> int:
    if args.samples <= 0:
        raise ValueError("--samples must be positive")
    if not args.unique_dual and not args.combinatorial:
        args.unique_dual = True
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    assignment = _assignment_from_args(args)
    i7_constraints = PanelConstraints(
        length=args.i7_length,
        min_hamming_distance=args.min_i7_distance,
        min_levenshtein_distance=max(1, args.min_i7_distance),
        min_sequence_levenshtein_distance=max(1, args.min_i7_distance),
        forbidden_motifs=[],
    )
    i5_constraints = PanelConstraints(
        length=args.i5_length,
        min_hamming_distance=args.min_i5_distance,
        min_levenshtein_distance=max(1, args.min_i5_distance),
        min_sequence_levenshtein_distance=max(1, args.min_i5_distance),
        forbidden_motifs=[],
    )
    i7_options = DesignOptions(
        panel_name="i7",
        count=args.i7_count,
        length=args.i7_length,
        seed=args.seed,
        restarts=5,
        candidate_pool_size=max(1500, args.i7_count * 250),
    )
    i5_options = DesignOptions(
        panel_name="i5",
        count=args.i5_count,
        length=args.i5_length,
        seed=args.seed + 101,
        restarts=5,
        candidate_pool_size=max(1500, args.i5_count * 250),
    )
    i7_records, _i7_trace = design_panel(i7_options, assignment, i7_constraints)
    i5_records, _i5_trace = design_panel(i5_options, assignment, i5_constraints)
    if len(i7_records) < args.i7_count or len(i5_records) < args.i5_count:
        raise ValueError("could not design enough i7/i5 barcodes")
    pair_rows: list[dict[str, object]] = []
    for idx in range(args.samples):
        i7 = i7_records[idx % len(i7_records)]
        i5 = i5_records[idx % len(i5_records)]
        pair_rows.append(
            {
                "sample_index": f"S{idx + 1:03d}",
                "i7_id": i7.barcode_id,
                "i7_seq": i7.sequence,
                "i5_id": i5.barcode_id,
                "i5_seq": i5.sequence,
                "concatenated_pair": f"{i7.sequence}+{i5.sequence}",
                "pair_min_distance": "",
                "i7_min_distance": "",
                "i5_min_distance": "",
                "index_hop_detectable": "true" if args.unique_dual else "false",
                "status": "pass",
            }
        )
    pair_sequences = [PanelRecord(str(row["sample_index"]), str(row["i7_seq"]) + str(row["i5_seq"])) for row in pair_rows]
    pair_min = _minimum_pairwise_distance(pair_sequences, "hamming")
    i7_min = _minimum_pairwise_distance(i7_records[: args.i7_count], "hamming")
    i5_min = _minimum_pairwise_distance(i5_records[: args.i5_count], "hamming")
    for row in pair_rows:
        row["pair_min_distance"] = pair_min if pair_min is not None else ""
        row["i7_min_distance"] = i7_min if i7_min is not None else ""
        row["i5_min_distance"] = i5_min if i5_min is not None else ""
        if isinstance(pair_min, int) and pair_min < args.min_pair_distance:
            row["status"] = "warn"
    _write_tsv(
        out_dir / "dual_barcodes.tsv",
        [
            "sample_index",
            "i7_id",
            "i7_seq",
            "i5_id",
            "i5_seq",
            "concatenated_pair",
            "pair_min_distance",
            "i7_min_distance",
            "i5_min_distance",
            "index_hop_detectable",
            "status",
        ],
        pair_rows,
    )
    _write_tsv(out_dir / "i7_barcodes.tsv", ["barcode_id", "sequence"], [{"barcode_id": r.barcode_id, "sequence": r.sequence} for r in i7_records[: args.i7_count]])
    _write_tsv(out_dir / "i5_barcodes.tsv", ["barcode_id", "sequence"], [{"barcode_id": r.barcode_id, "sequence": r.sequence} for r in i5_records[: args.i5_count]])
    check_panel(i7_records[: args.i7_count], assignment, i7_constraints, out_dir / "panel_check_i7", source_panel=out_dir / "i7_barcodes.tsv")
    check_panel(i5_records[: args.i5_count], assignment, i5_constraints, out_dir / "panel_check_i5", source_panel=out_dir / "i5_barcodes.tsv")
    export_dual_panel(pair_rows, assignment, out_dir)
    _write_json(
        out_dir / "design_report.json",
        {
            "schema_version": PANEL_SCHEMA_VERSION,
            "mode": "design-dual",
            "samples": args.samples,
            "unique_dual": bool(args.unique_dual),
            "combinatorial": bool(args.combinatorial),
            "i7_count": args.i7_count,
            "i5_count": args.i5_count,
            "i7_min_distance": i7_min,
            "i5_min_distance": i5_min,
            "pair_min_distance": pair_min,
            "index_hop_detectable": bool(args.unique_dual),
            "certified_dotmatch_command": _certified_command(out_dir / "dual_barcodes.tsv", assignment, args.i7_length + args.i5_length),
        },
    )
    _write_basic_html(
        out_dir / "report.html",
        "DotMatch Dual-Index Panel Report",
        [
            f"Samples: {args.samples}",
            f"i7 minimum Hamming distance: {i7_min}",
            f"i5 minimum Hamming distance: {i5_min}",
            f"Pair minimum Hamming distance: {pair_min}",
            "Index-hop detectable: yes" if args.unique_dual else "Index-hop detectable: combinatorial mode is higher risk",
        ],
    )
    return 0


def read_panel(path: str | Path) -> list[PanelRecord]:
    p = Path(path)
    delimiter = "," if p.suffix.lower() == ".csv" else "\t"
    lines = [line for line in p.read_text(encoding="utf-8").splitlines() if line.strip() and not line.startswith("#")]
    if not lines:
        raise ValueError(f"no barcodes found in {path}")
    first = lines[0].split(delimiter)
    normalized = [col.strip().lower() for col in first]
    has_header = bool(set(normalized) & {"barcode_id", "id", "sample_id", "sample_index"}) and bool(
        set(normalized) & {"sequence", "seq", "barcode_seq", "i7_seq", "index"}
    )
    records: list[PanelRecord] = []
    if has_header:
        reader = csv.DictReader(lines, delimiter=delimiter)
        assert reader.fieldnames is not None
        id_col = _first_present(reader.fieldnames, ["barcode_id", "id", "sample_id", "sample_index", "target_id"]) or reader.fieldnames[0]
        seq_col = _first_present(reader.fieldnames, ["sequence", "barcode_seq", "seq", "i7_seq", "index", "target_seq"])
        if seq_col is None:
            raise ValueError("panel file must contain a sequence column")
        for i, row in enumerate(reader):
            seq = str(row.get(seq_col, "") or "").strip().upper()
            if not seq:
                continue
            barcode_id = str(row.get(id_col, "") or f"BC{i + 1:03d}").strip()
            records.append(PanelRecord(barcode_id or f"BC{i + 1:03d}", seq))
    else:
        for i, line in enumerate(lines):
            cols = line.split(delimiter)
            if len(cols) == 1:
                records.append(PanelRecord(f"BC{i + 1:03d}", cols[0].strip().upper()))
            else:
                records.append(PanelRecord(cols[0].strip() or f"BC{i + 1:03d}", cols[1].strip().upper()))
    if not records:
        raise ValueError(f"no barcodes found in {path}")
    return records


def check_panel(
    records: Sequence[PanelRecord],
    assignment: AssignmentConfig,
    constraints: PanelConstraints,
    out_dir: Path,
    *,
    source_panel: Path | None = None,
) -> dict[str, object]:
    if assignment.k < 0:
        raise ValueError("--k must be non-negative")
    if assignment.k > MAX_EXACT_CERTIFICATE_RADIUS:
        raise ValueError(f"exact safety certificate supports k <= {MAX_EXACT_CERTIFICATE_RADIUS}; refusing to certify k={assignment.k}")
    out_dir.mkdir(parents=True, exist_ok=True)
    source_panel = source_panel or Path("barcodes.tsv")
    lengths = sorted({len(record.sequence) for record in records})
    pair_rows, nearest = _pairwise_rows(records, constraints)
    duplicate_pairs = [row for row in pair_rows if row["hamming_distance"] == 0 and row["levenshtein_distance"] == 0]
    prefix_rows = _prefix_suffix_rows(records, assignment, prefix=True)
    suffix_rows = _prefix_suffix_rows(records, assignment, prefix=False)
    rc_rows = _reverse_complement_rows(records, assignment, constraints)
    content_rows = _content_rows(records, constraints, nearest)
    cycle_rows, cycle_warning_count, cycle_score = _cycle_balance_rows(records, constraints)
    context_rows = _context_rows(records, constraints)
    flanked_rows = _flanked_rows(records, constraints)

    configured_risk = _assignment_risk_rows(records, assignment, constraints)
    ambiguous_rows = [row for row in configured_risk if row["risk_type"] == "ambiguous"]
    silent_rows = [row for row in configured_risk if row["risk_type"] == "silent_assignment"]
    target_safety_rows = _target_safety_rows(records, content_rows, configured_risk)

    safe_k0 = len(duplicate_pairs) == 0
    safe_k1_hamming = _safe_for_radius(records, "hamming", 1, constraints)
    safe_k1_levenshtein = _safe_for_radius(records, "levenshtein", 1, constraints)
    safe_k2_hamming = _safe_by_distance_rule(records, "hamming", 2, constraints)
    safe_k2_levenshtein = _safe_by_distance_rule(records, "levenshtein", 2, constraints)
    configured_safe = not ambiguous_rows and not silent_rows
    has_prefix_fail = bool(prefix_rows) and constraints.avoid_prefix_collisions
    has_suffix_fail = bool(suffix_rows) and constraints.avoid_suffix_collisions
    has_rc_fail = bool(rc_rows) and assignment.reverse_complement_mode == "fail"
    content_fail = any(row["status"] == "fail" for row in content_rows)
    status = "pass"
    if duplicate_pairs or not configured_safe or has_prefix_fail or has_suffix_fail or has_rc_fail or content_fail:
        status = "fail"
    elif rc_rows or cycle_warning_count or context_rows or _narrow_margin(records, assignment):
        status = "warn"

    min_hamming = _minimum_pairwise_distance(records, "hamming")
    min_lev = _minimum_pairwise_distance(records, "levenshtein")
    min_seqlev = _minimum_pairwise_distance(records, "sequence-levenshtein", constraints)
    certified_command = _certified_command(source_panel, assignment, lengths[0] if len(lengths) == 1 else "auto")
    for rows in [pair_rows, prefix_rows, suffix_rows, rc_rows, content_rows, cycle_rows, context_rows, flanked_rows, ambiguous_rows, configured_risk, target_safety_rows]:
        _stamp_command(rows, certified_command)
    summary = {
        "schema_version": PANEL_SCHEMA_VERSION,
        "status": status,
        "panel_grade": _panel_grade(status, records, assignment, constraints, min_hamming, min_lev, ambiguous_rows, silent_rows, rc_rows, prefix_rows, suffix_rows),
        "n_barcodes": len(records),
        "lengths": lengths,
        "assignment_metric": _normalize_metric(assignment.metric),
        "configured_assignment_k": assignment.k,
        "exact_error_sphere_radius": assignment.k,
        "maximum_exact_certificate_radius": MAX_EXACT_CERTIFICATE_RADIUS,
        "ambiguity_policy": assignment.ambiguity_policy,
        "safe_for_k0": safe_k0,
        "safe_for_k1_hamming": safe_k1_hamming,
        "safe_for_k1_levenshtein": safe_k1_levenshtein,
        "safe_for_k2_hamming": safe_k2_hamming,
        "safe_for_k2_levenshtein": safe_k2_levenshtein,
        "minimum_hamming_distance": min_hamming,
        "minimum_levenshtein_distance": min_lev,
        "minimum_sequence_levenshtein_distance": min_seqlev,
        "guaranteed_substitution_correction": _guaranteed_correction(min_hamming),
        "guaranteed_substitution_detection": _guaranteed_detection(min_hamming),
        "ambiguous_error_spheres": len(ambiguous_rows),
        "unsafe_k1_variants": len(configured_risk) if assignment.k == 1 else len(ambiguous_rows) + len(silent_rows),
        "unsafe_configured_radius_variants": len(configured_risk),
        "silent_assignment_risk": len(silent_rows),
        "duplicate_pairs": len(duplicate_pairs),
        "collision_pairs": sum(1 for row in pair_rows if row["status"] != "pass"),
        "prefix_collisions": len(prefix_rows),
        "suffix_collisions": len(suffix_rows),
        "reverse_complement_warnings": len(rc_rows),
        "cycle_balance_warnings": cycle_warning_count,
        "cycle_balance_score": cycle_score,
        "forbidden_motif_hits": sum(1 for row in content_rows if row["forbidden_motifs"]),
        "maximum_homopolymer": max((int(row["homopolymer_max"]) for row in content_rows), default=0),
        "certified_dotmatch_command": certified_command,
        "scope": "barcode design for known-target assignment; not general genome alignment, not UMI entropy generation, not basecalling",
    }

    _write_json(out_dir / "panel_summary.json", summary)
    _write_tsv(out_dir / "panel_summary.tsv", list(summary.keys()), [summary])
    _write_tsv(out_dir / "target_safety.tsv", _target_safety_columns(), target_safety_rows)
    _write_tsv(out_dir / "collision_pairs.tsv", _pairwise_columns(), pair_rows)
    _write_tsv(out_dir / "ambiguous_error_spheres.tsv", _risk_columns(), ambiguous_rows)
    _write_tsv(out_dir / "unsafe_assignment_variants.tsv", _risk_columns(), configured_risk)
    _write_tsv(out_dir / "prefix_collisions.tsv", _prefix_columns(), prefix_rows)
    _write_tsv(out_dir / "suffix_collisions.tsv", _prefix_columns(), suffix_rows)
    _write_tsv(out_dir / "reverse_complement_warnings.tsv", _rc_columns(), rc_rows)
    _write_tsv(out_dir / "homopolymer_gc.tsv", _content_columns(), content_rows)
    _write_tsv(out_dir / "cycle_balance.tsv", _cycle_columns(), cycle_rows)
    _write_tsv(out_dir / "context_risk.tsv", _context_columns(), context_rows)
    _write_tsv(out_dir / "flanked_sequences.tsv", _flanked_columns(), flanked_rows)
    _write_panel_report(out_dir / "panel_report.html", summary, records, pair_rows, ambiguous_rows, rc_rows, prefix_rows, suffix_rows)
    _write_panel_report(out_dir / "report.html", summary, records, pair_rows, ambiguous_rows, rc_rows, prefix_rows, suffix_rows)
    return summary


def design_panel(
    options: DesignOptions,
    assignment: AssignmentConfig,
    constraints: PanelConstraints,
) -> tuple[list[PanelRecord], list[dict[str, object]]]:
    if options.count <= 0:
        raise ValueError("--n must be positive")
    if options.length <= 0:
        raise ValueError("--length must be positive")
    constraints.length = options.length
    rng = random.Random(options.seed)
    candidates = generate_candidate_pool(options, constraints, rng)
    if not candidates:
        raise ValueError("no candidates survived hard filters")
    if options.engine == "graph":
        return select_panel_graph(candidates, options, constraints)
    if options.engine == "evolve":
        return select_panel_evolve(candidates, options, constraints)
    return select_panel_from_candidates(candidates, options, constraints)


def generate_candidate_pool(options: DesignOptions, constraints: PanelConstraints, rng: random.Random) -> list[PanelRecord]:
    seen: set[str] = set()
    candidates: list[PanelRecord] = []
    max_attempts = max(options.candidate_pool_size * 8, options.count * 500)
    attempts = 0
    while attempts < max_attempts and len(candidates) < options.candidate_pool_size:
        attempts += 1
        seq = "".join(rng.choice(DNA) for _ in range(options.length))
        if seq in seen:
            continue
        seen.add(seq)
        if candidate_filter_reason(seq, constraints) != "pass":
            continue
        candidates.append(PanelRecord(f"CAND{len(candidates) + 1:06d}", seq))
    candidates.sort(key=lambda record: (-individual_candidate_score(record.sequence, constraints), record.sequence))
    return candidates


def select_panel_from_candidates(
    candidates: Sequence[PanelRecord],
    options: DesignOptions,
    constraints: PanelConstraints,
) -> tuple[list[PanelRecord], list[dict[str, object]]]:
    rng = random.Random(options.seed)
    best: list[PanelRecord] = []
    best_score = -float("inf")
    best_trace: list[dict[str, object]] = []
    start = time.monotonic()
    working = list(candidates)
    working.sort(key=lambda record: (-individual_candidate_score(record.sequence, constraints), record.sequence))
    search_limit = min(len(working), max(options.count * 1000, 1200))
    working = working[:search_limit]
    restarts = max(1, options.restarts)
    for restart in range(restarts):
        if options.max_runtime and time.monotonic() - start > options.max_runtime:
            break
        selected: list[PanelRecord] = []
        trace: list[dict[str, object]] = []
        restart_pool = working[:]
        rng.shuffle(restart_pool)
        restart_pool.sort(key=lambda record: -(individual_candidate_score(record.sequence, constraints) + rng.random() * 0.0001))
        while len(selected) < options.count:
            best_candidate: PanelRecord | None = None
            best_candidate_score = -float("inf")
            checked = 0
            for candidate in restart_pool:
                if candidate in selected:
                    continue
                if not all(pair_compatible(candidate.sequence, item.sequence, constraints) for item in selected):
                    continue
                score = incremental_candidate_score(candidate, selected, constraints, options.objective)
                checked += 1
                if score > best_candidate_score:
                    best_candidate = candidate
                    best_candidate_score = score
                if checked > max(500, options.count * 40):
                    break
            if best_candidate is None:
                break
            selected.append(best_candidate)
            trace.append(
                {
                    "restart": restart,
                    "step": len(selected),
                    "candidate_id": best_candidate.barcode_id,
                    "sequence": best_candidate.sequence,
                    "score": f"{best_candidate_score:.6f}",
                    "selected_count": len(selected),
                    "reason": "selected",
                }
            )
        score = panel_objective(selected, constraints, options.objective) + len(selected) * 10000
        if len(selected) > len(best) or (len(selected) == len(best) and score > best_score):
            best = selected
            best_score = score
            best_trace = trace
        if len(best) >= options.count and restart >= 1:
            # Once a full panel is found, a few restarts are enough for deterministic v1 output.
            continue
    return _renumber_records(best[: options.count]), best_trace


def select_panel_graph(
    candidates: Sequence[PanelRecord],
    options: DesignOptions,
    constraints: PanelConstraints,
) -> tuple[list[PanelRecord], list[dict[str, object]]]:
    # Greedy clique over a compatibility graph. This is intentionally exact about edge
    # compatibility and heuristic about clique search, which keeps v1 predictable.
    top = list(candidates)[: min(len(candidates), max(800, options.count * 250))]
    top.sort(key=lambda record: (-individual_candidate_score(record.sequence, constraints), record.sequence))
    selected: list[PanelRecord] = []
    trace: list[dict[str, object]] = []
    start = time.monotonic()
    while len(selected) < options.count:
        if options.max_runtime and time.monotonic() - start > options.max_runtime:
            break
        compatible = [record for record in top if record not in selected and all(pair_compatible(record.sequence, other.sequence, constraints) for other in selected)]
        if not compatible:
            break
        compatible.sort(key=lambda record: (-incremental_candidate_score(record, selected, constraints, options.objective), record.sequence))
        chosen = compatible[0]
        selected.append(chosen)
        trace.append(
            {
                "restart": 0,
                "step": len(selected),
                "candidate_id": chosen.barcode_id,
                "sequence": chosen.sequence,
                "score": f"{panel_objective(selected, constraints, options.objective):.6f}",
                "selected_count": len(selected),
                "reason": "graph_clique_selected",
            }
        )
    return _renumber_records(selected[: options.count]), trace


def select_panel_evolve(
    candidates: Sequence[PanelRecord],
    options: DesignOptions,
    constraints: PanelConstraints,
) -> tuple[list[PanelRecord], list[dict[str, object]]]:
    base_options = DesignOptions(**{**options.__dict__, "engine": "greedy", "restarts": max(3, min(options.restarts, 10))})
    best, trace = select_panel_from_candidates(candidates, base_options, constraints)
    rng = random.Random(options.seed + 13)
    if len(best) < options.count:
        return best, trace
    best_score = panel_objective(best, constraints, options.objective)
    candidate_list = list(candidates)
    max_iterations = max(1, options.iterations)
    for iteration in range(max_iterations):
        trial = best[:]
        if not trial:
            break
        replace_index = rng.randrange(len(trial))
        rng.shuffle(candidate_list)
        for candidate in candidate_list[: max(1000, options.count * 200)]:
            if candidate.sequence in {record.sequence for record in trial}:
                continue
            proposal = trial[:]
            proposal[replace_index] = candidate
            if all(pair_compatible(a.sequence, b.sequence, constraints) for a, b in combinations(proposal, 2)):
                score = panel_objective(proposal, constraints, options.objective)
                if score > best_score or rng.random() < options.mutation_rate * 0.01:
                    best = proposal
                    best_score = score
                    trace.append(
                        {
                            "restart": "evolve",
                            "step": iteration + 1,
                            "candidate_id": candidate.barcode_id,
                            "sequence": candidate.sequence,
                            "score": f"{score:.6f}",
                            "selected_count": len(best),
                            "reason": "evolutionary_replacement",
                        }
                    )
                break
    return _renumber_records(best[: options.count]), trace


def select_panel_from_candidates_for_export(candidates: Sequence[PanelRecord], n: int, constraints: PanelConstraints, seed: int = 1) -> list[PanelRecord]:
    options = DesignOptions(count=n, length=constraints.length or len(candidates[0].sequence), seed=seed, candidate_pool_size=len(candidates))
    selected, _trace = select_panel_from_candidates(candidates, options, constraints)
    return selected


def simulate_panel(
    records: Sequence[PanelRecord],
    assignment: AssignmentConfig,
    *,
    reads: int,
    seed: int,
    substitution_rate: float,
    insertion_rate: float,
    deletion_rate: float,
    quality_model: str,
    out_dir: Path,
    source_panel: Path,
) -> dict[str, object]:
    if reads <= 0:
        raise ValueError("--reads must be positive")
    rng = random.Random(seed)
    assignment_rows: list[dict[str, object]] = []
    false_rows: list[dict[str, object]] = []
    ambiguous_examples: list[dict[str, object]] = []
    none_examples: list[dict[str, object]] = []
    confusion: dict[str, Counter[str]] = {record.barcode_id: Counter() for record in records}
    counts = Counter()
    for idx in range(reads):
        source_index = rng.randrange(len(records))
        source = records[source_index]
        observed = mutate_sequence(source.sequence, rng, substitution_rate, insertion_rate, deletion_rate)
        outcome = assign_query(observed, records, assignment.metric, assignment.k)
        assigned_id = records[outcome.target_index].barcode_id if outcome.status == "unique" and outcome.target_index >= 0 else ""
        assigned_sequence = records[outcome.target_index].sequence if outcome.status == "unique" and outcome.target_index >= 0 else ""
        false_assignment = outcome.status == "unique" and outcome.target_index != source_index
        counts[outcome.status] += 1
        if false_assignment:
            counts["false_assignment"] += 1
        confusion[source.barcode_id][assigned_id or outcome.status] += 1
        row = {
            "read_id": f"sim_{idx + 1:08d}",
            "true_id": source.barcode_id,
            "true_sequence": source.sequence,
            "observed_sequence": observed,
            "assigned_id": assigned_id,
            "assigned_sequence": assigned_sequence,
            "distance": outcome.best_distance,
            "status": outcome.status,
            "false_assignment": "true" if false_assignment else "false",
            "certified_command": _certified_command(source_panel, assignment, _single_length(records) or "auto"),
        }
        assignment_rows.append(row)
        if false_assignment:
            false_rows.append(row)
        if outcome.status == "ambiguous" and len(ambiguous_examples) < 100:
            ambiguous_examples.append(row)
        if outcome.status == "none" and len(none_examples) < 100:
            none_examples.append(row)
    risk_rows = _assignment_risk_rows(records, assignment, PanelConstraints(length=_single_length(records)))
    per_barcode_rows: list[dict[str, object]] = []
    worst_recall = 1.0
    worst_pair = ""
    worst_pair_count = 0
    for record in records:
        source_total = sum(confusion[record.barcode_id].values())
        true_positive = confusion[record.barcode_id].get(record.barcode_id, 0)
        assigned_total = sum(confusion[source.barcode_id].get(record.barcode_id, 0) for source in records)
        recall = true_positive / source_total if source_total else 0.0
        precision = true_positive / assigned_total if assigned_total else 0.0
        worst_recall = min(worst_recall, recall)
        for assigned_key, count in confusion[record.barcode_id].items():
            if assigned_key != record.barcode_id and count > worst_pair_count:
                worst_pair_count = count
                worst_pair = f"{record.barcode_id}->{assigned_key}"
        per_barcode_rows.append(
            {
                "barcode_id": record.barcode_id,
                "true_reads": source_total,
                "unique_correct": true_positive,
                "recall": f"{recall:.8f}",
                "precision": f"{precision:.8f}",
                "most_common_outcome": confusion[record.barcode_id].most_common(1)[0][0] if source_total else "",
                "certified_command": _certified_command(source_panel, assignment, _single_length(records) or "auto"),
            }
        )
    false_rate = counts["false_assignment"] / reads
    summary = {
        "schema_version": PANEL_SCHEMA_VERSION,
        "total_reads": reads,
        "seed": seed,
        "quality_model": quality_model,
        "assignment_metric": _normalize_metric(assignment.metric),
        "k": assignment.k,
        "unique": counts["unique"],
        "ambiguous": counts["ambiguous"],
        "none": counts["none"],
        "invalid": counts["invalid"],
        "false_assignment": counts["false_assignment"],
        "unique_rate": counts["unique"] / reads,
        "ambiguous_rate": counts["ambiguous"] / reads,
        "none_rate": counts["none"] / reads,
        "invalid_rate": counts["invalid"] / reads,
        "false_assignment_rate": false_rate,
        "false_assignment_upper_bound": false_rate + 1.96 * math.sqrt(false_rate * (1 - false_rate) / reads) if reads else 0.0,
        "worst_barcode_recall": worst_recall,
        "worst_pair_confusion": worst_pair,
        "ambiguous_variant_count": sum(1 for row in risk_rows if row["risk_type"] == "ambiguous"),
        "unsafe_variant_count": len(risk_rows),
        "certified_dotmatch_command": _certified_command(source_panel, assignment, _single_length(records) or "auto"),
    }
    _write_tsv(out_dir / "simulated_assignments.tsv", _simulation_columns(), assignment_rows)
    _write_json(out_dir / "simulation_summary.json", summary)
    _write_tsv(out_dir / "false_assignment.tsv", _simulation_columns(), false_rows)
    _write_tsv(out_dir / "ambiguous_examples.tsv", _simulation_columns(), ambiguous_examples)
    _write_tsv(out_dir / "none_examples.tsv", _simulation_columns(), none_examples)
    _write_tsv(out_dir / "per_barcode_confusion.tsv", _confusion_columns(), per_barcode_rows)
    _write_basic_html(
        out_dir / "simulation_report.html",
        "DotMatch Panel Simulation Report",
        [
            f"Unique rate: {summary['unique_rate']:.4f}",
            f"Ambiguous rate: {summary['ambiguous_rate']:.4f}",
            f"None rate: {summary['none_rate']:.4f}",
            f"False assignment rate: {summary['false_assignment_rate']:.6f}",
            f"Command: {summary['certified_dotmatch_command']}",
        ],
    )
    return summary


def write_plate_layout(records: Sequence[PanelRecord], assignment: AssignmentConfig, plate: int, out: Path) -> None:
    rows, cols = _plate_shape(plate)
    wells = [f"{chr(ord('A') + r)}{c + 1}" for r in range(rows) for c in range(cols)]
    if len(records) > len(wells):
        raise ValueError(f"panel has {len(records)} barcodes but plate {plate} has {len(wells)} wells")
    ordered = _plate_order(records)
    placed: dict[str, PanelRecord] = {}
    output_rows: list[dict[str, object]] = []
    for well, record in zip(wells, ordered):
        placed[well] = record
        row_letter = "".join(ch for ch in well if ch.isalpha())
        column = int("".join(ch for ch in well if ch.isdigit()))
        output_rows.append(
            {
                "well": well,
                "row": row_letter,
                "column": column,
                "barcode_id": record.barcode_id,
                "sequence": record.sequence,
                "nearest_neighbor_distance": _nearest_distance(record, [item for item in records if item != record], "hamming"),
                "certified_command": _certified_command(out, assignment, _single_length(records) or "auto"),
            }
        )
    _write_tsv(out, ["well", "row", "column", "barcode_id", "sequence", "nearest_neighbor_distance", "certified_command"], output_rows)
    neighbor_rows = _neighbor_rows(output_rows, records)
    _write_tsv(out.parent / "neighbor_distance.tsv", ["well", "neighbor_well", "barcode_id", "neighbor_barcode_id", "hamming_distance"], neighbor_rows)
    _write_picklist(out.parent / "lab_picklist.csv", output_rows)
    _write_plate_svg(out.parent / f"{out.stem}.svg", output_rows, rows, cols)


def export_panel(records: Sequence[PanelRecord], assignment: AssignmentConfig, out_dir: Path, *, fmt: str, source_panel: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    if fmt == "illumina-samplesheet":
        _write_illumina_samplesheet(out_dir / "SampleSheet.csv", records)
    else:
        _write_tsv(out_dir / "panel_export.tsv", ["sample_id", "index", "barcode_id", "certified_command"], [
            {
                "sample_id": record.barcode_id,
                "index": record.sequence,
                "barcode_id": record.barcode_id,
                "certified_command": _certified_command(source_panel, assignment, _single_length(records) or "auto"),
            }
            for record in records
        ])
    _write_tsv(out_dir / "barcodes_for_demux.tsv", ["barcode_id", "sequence", "certified_command"], [
        {
            "barcode_id": record.barcode_id,
            "sequence": record.sequence,
            "certified_command": _certified_command(source_panel, assignment, _single_length(records) or "auto"),
        }
        for record in records
    ])
    _write_lab_readme(out_dir / "README_FOR_LAB.md", records, assignment, source_panel)


def export_dual_panel(pair_rows: Sequence[dict[str, object]], assignment: AssignmentConfig, out_dir: Path) -> None:
    templates = out_dir / "sample_sheet_templates"
    templates.mkdir(parents=True, exist_ok=True)
    lines = ["[Header]", "IEMFileVersion,5", "Workflow,GenerateFASTQ", "", "[Data]", "Sample_ID,index,index2"]
    for row in pair_rows:
        lines.append(f"{row['sample_index']},{row['i7_seq']},{row['i5_seq']}")
    (templates / "SampleSheet.csv").write_text("\n".join(lines) + "\n", encoding="utf-8")
    _write_lab_readme(out_dir / "README_FOR_LAB.md", [], assignment, out_dir / "dual_barcodes.tsv")


def compare_panels(
    old_records: Sequence[PanelRecord],
    new_records: Sequence[PanelRecord],
    assignment: AssignmentConfig,
    out_dir: Path,
    old_panel: Path,
    new_panel: Path,
) -> None:
    old_sequences = {record.sequence for record in old_records}
    new_sequences = {record.sequence for record in new_records}
    old_ids = {record.barcode_id for record in old_records}
    new_ids = {record.barcode_id for record in new_records}
    summary = {
        "schema_version": PANEL_SCHEMA_VERSION,
        "old_count": len(old_records),
        "new_count": len(new_records),
        "shared_sequences": len(old_sequences & new_sequences),
        "added_sequences": len(new_sequences - old_sequences),
        "removed_sequences": len(old_sequences - new_sequences),
        "shared_ids": len(old_ids & new_ids),
        "old_min_hamming_distance": _minimum_pairwise_distance(old_records, "hamming"),
        "new_min_hamming_distance": _minimum_pairwise_distance(new_records, "hamming"),
        "old_panel": str(old_panel),
        "new_panel": str(new_panel),
        "certified_dotmatch_command": _certified_command(new_panel, assignment, _single_length(new_records) or "auto"),
    }
    _write_json(out_dir / "panel_compare.json", summary)
    rows = []
    for seq in sorted(new_sequences - old_sequences):
        rows.append({"change": "added", "sequence": seq})
    for seq in sorted(old_sequences - new_sequences):
        rows.append({"change": "removed", "sequence": seq})
    _write_tsv(out_dir / "panel_compare.tsv", ["change", "sequence"], rows)


def _write_design_outputs(
    out_dir: Path,
    records: Sequence[PanelRecord],
    trace: Sequence[dict[str, object]],
    options: DesignOptions,
    assignment: AssignmentConfig,
    constraints: PanelConstraints,
    *,
    mode: str,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    temp_check_dir = out_dir / "panel_check"
    certificate = check_panel(records, assignment, constraints, temp_check_dir, source_panel=Path("barcodes.tsv"))
    write_barcode_table(out_dir / "barcodes.tsv", records, constraints, certificate)
    # Re-run after barcodes.tsv exists so report commands point at the final file.
    certificate = check_panel(records, assignment, constraints, temp_check_dir, source_panel=Path("barcodes.tsv"))
    write_barcode_table(out_dir / "barcodes.tsv", records, constraints, certificate)
    if certificate["status"] == "fail":
        raise ValueError("designed panel failed its DotMatch safety certificate; adjust constraints, seed, or candidate pool")
    _copy_if_exists(temp_check_dir / "target_safety.tsv", out_dir / "assignment_safety.tsv")
    _copy_if_exists(temp_check_dir / "collision_pairs.tsv", out_dir / "collision_pairs.tsv")
    _copy_if_exists(temp_check_dir / "ambiguous_error_spheres.tsv", out_dir / "ambiguous_error_spheres.tsv")
    _copy_if_exists(temp_check_dir / "panel_report.html", out_dir / "report.html")
    write_plate_layout(records, assignment, options.plate_format, out_dir / "plate_layout.tsv")
    sample_dir = out_dir / "sample_sheet_templates"
    export_panel(records, assignment, sample_dir, fmt="illumina-samplesheet", source_panel=out_dir / "barcodes.tsv")
    _write_lab_readme(out_dir / "README_FOR_LAB.md", records, assignment, out_dir / "barcodes.tsv")
    _write_tsv(out_dir / "design_trace.tsv", ["restart", "step", "candidate_id", "sequence", "score", "selected_count", "reason"], trace)
    design_report = {
        "schema_version": PANEL_SCHEMA_VERSION,
        "mode": mode,
        "panel_name": options.panel_name,
        "engine": options.engine,
        "objective": options.objective,
        "weights": OBJECTIVE_WEIGHTS[options.objective],
        "seed": options.seed,
        "requested_count": options.count,
        "selected_count": len(records),
        "length": options.length,
        "candidate_pool_size": options.candidate_pool_size,
        "restarts": options.restarts,
        "certificate": certificate,
        "certified_dotmatch_command": certificate["certified_dotmatch_command"],
    }
    _write_json(out_dir / "design_report.json", design_report)
    if options.experimental_signal or options.objective == "ont-signal-separable":
        _write_signal_outputs(out_dir, records, options, certificate)


def write_barcode_table(path: Path, records: Sequence[PanelRecord], constraints: PanelConstraints, certificate: dict[str, object]) -> None:
    nearest: dict[str, dict[str, object]] = {}
    for record in records:
        others = [item for item in records if item != record]
        nearest[record.barcode_id] = {
            "hamming": _nearest_distance(record, others, "hamming"),
            "levenshtein": _nearest_distance(record, others, "levenshtein"),
            "seqlev": _nearest_distance(record, others, "sequence-levenshtein", constraints),
        }
    rc_sequences = {reverse_complement(record.sequence) for record in records}
    rows = []
    for record in records:
        content = sequence_content(record.sequence, constraints)
        rows.append(
            {
                "barcode_id": record.barcode_id,
                "sequence": record.sequence,
                "length": len(record.sequence),
                "gc": f"{content['gc']:.2f}",
                "homopolymer_max": content["homopolymer_max"],
                "min_hamming_neighbor": nearest[record.barcode_id]["hamming"],
                "min_lev_neighbor": nearest[record.barcode_id]["levenshtein"],
                "min_seqlev_neighbor": nearest[record.barcode_id]["seqlev"],
                "rc_collision": "true" if reverse_complement(record.sequence) in rc_sequences - {record.sequence} else "false",
                "self_complement_score": f"{self_complement_score(record.sequence):.2f}",
                "status": "pass" if content["status"] == "pass" else content["status"],
                "warnings": ";".join(content["warnings"]),
                "certified_command": certificate["certified_dotmatch_command"],
            }
        )
    _write_tsv(
        path,
        [
            "barcode_id",
            "sequence",
            "length",
            "gc",
            "homopolymer_max",
            "min_hamming_neighbor",
            "min_lev_neighbor",
            "min_seqlev_neighbor",
            "rc_collision",
            "self_complement_score",
            "status",
            "warnings",
            "certified_command",
        ],
        rows,
    )


def _assignment_from_args(args: argparse.Namespace) -> AssignmentConfig:
    if getattr(args, "k", 0) < 0:
        raise ValueError("--k must be non-negative")
    return AssignmentConfig(
        metric=_normalize_metric(args.metric),
        k=args.k,
        reverse_complement_mode=args.reverse_complement_mode,
    )


def _constraints_from_args(args: argparse.Namespace, *, length: int | None = None) -> PanelConstraints:
    constraints = PanelConstraints(length=length)
    if getattr(args, "min_hamming_distance", None) is not None:
        constraints.min_hamming_distance = args.min_hamming_distance
    if getattr(args, "min_levenshtein_distance", None) is not None:
        constraints.min_levenshtein_distance = args.min_levenshtein_distance
    if getattr(args, "min_sequence_levenshtein_distance", None) is not None:
        constraints.min_sequence_levenshtein_distance = args.min_sequence_levenshtein_distance
    if getattr(args, "gc_min", None) is not None:
        constraints.gc_min = args.gc_min
    if getattr(args, "gc_max", None) is not None:
        constraints.gc_max = args.gc_max
    if getattr(args, "max_homopolymer", None) is not None:
        constraints.max_homopolymer = args.max_homopolymer
    if getattr(args, "max_dinucleotide_repeat", None) is not None:
        constraints.max_dinucleotide_repeat = args.max_dinucleotide_repeat
    if getattr(args, "avoid_reverse_complements", None) is not None:
        constraints.avoid_reverse_complements = bool(args.avoid_reverse_complements)
    if getattr(args, "allow_ambiguous_literals", False):
        constraints.allow_ambiguous_literals = True
        constraints.avoid_ambiguous_bases = False
    if getattr(args, "forbidden_motifs", None):
        constraints.forbidden_motifs = [motif.upper() for motif in args.forbidden_motifs]
    if hasattr(args, "left_flank"):
        constraints.left_flank = str(args.left_flank or "").upper()
    if hasattr(args, "right_flank"):
        constraints.right_flank = str(args.right_flank or "").upper()
    if getattr(args, "context_window", None) is not None:
        constraints.context_window = args.context_window
    return constraints


def _design_from_args(args: argparse.Namespace) -> tuple[DesignOptions, AssignmentConfig, PanelConstraints]:
    spec_data = _load_design_spec(Path(args.spec)) if args.spec else {}
    panel_data = _as_dict(spec_data.get("panel"))
    assignment_data = _as_dict(spec_data.get("assignment"))
    constraints_data = _as_dict(spec_data.get("constraints"))
    cycle_data = _as_dict(spec_data.get("cycle_balance"))
    plate_data = _as_dict(spec_data.get("plate_layout"))
    simulation_data = _as_dict(spec_data.get("simulation"))
    merged: dict[str, object] = {}
    if args.preset:
        merged.update(PRESETS[args.preset])
    merged.update(
        {
            "panel_name": panel_data.get("name", merged.get("panel_name", "dotmatch_panel")),
            "count": panel_data.get("count", panel_data.get("n", merged.get("count", 96))),
            "length": panel_data.get("length", merged.get("length", 16)),
            "seed": panel_data.get("seed", merged.get("seed", 1)),
            "alphabet": panel_data.get("alphabet", DNA),
            "metric": assignment_data.get("metric", "hamming"),
            "k": assignment_data.get("k", 1),
            "reverse_complement_mode": assignment_data.get("reverse_complement_mode", "warn"),
            "plate_format": plate_data.get("format", 96),
        }
    )
    merged.update(constraints_data)
    if cycle_data:
        merged["cycle_balance_enabled"] = cycle_data.get("enabled", True)
        merged["cycle_min_base_fraction"] = cycle_data.get("min_base_fraction_per_cycle", 0.15)
        merged["cycle_max_base_fraction"] = cycle_data.get("max_base_fraction_per_cycle", 0.40)
    if simulation_data:
        merged["simulation"] = simulation_data
    if args.count is not None:
        merged["count"] = args.count
    if args.length is not None:
        merged["length"] = args.length
    if args.seed is not None:
        merged["seed"] = args.seed
    if args.engine is not None:
        merged["engine"] = args.engine
    if args.restarts is not None:
        merged["restarts"] = args.restarts
    if args.candidate_pool_size is not None:
        merged["candidate_pool_size"] = args.candidate_pool_size
    if args.max_runtime is not None:
        merged["max_runtime"] = args.max_runtime
    if args.target_count is not None:
        merged["target_count"] = args.target_count
    if args.population is not None:
        merged["population"] = args.population
    if args.iterations is not None:
        merged["iterations"] = args.iterations
    if args.mutation_rate is not None:
        merged["mutation_rate"] = args.mutation_rate
    if args.elite_fraction is not None:
        merged["elite_fraction"] = args.elite_fraction
    if args.objective is not None:
        merged["objective"] = args.objective
    if getattr(args, "metric", None) is not None:
        merged["metric"] = args.metric
    if getattr(args, "k", None) is not None:
        merged["k"] = args.k
    if getattr(args, "reverse_complement_mode", None) is not None:
        merged["reverse_complement_mode"] = args.reverse_complement_mode

    options = DesignOptions(
        panel_name=str(merged.get("panel_name", "dotmatch_panel")),
        count=int(merged.get("count", 96)),
        length=int(merged.get("length", 16)),
        seed=int(merged.get("seed", 1)),
        engine=str(merged.get("engine", "greedy")),
        restarts=int(merged.get("restarts", 100)),
        candidate_pool_size=int(merged.get("candidate_pool_size", 100000)),
        objective=str(merged.get("objective", "strict-demux")),
        max_runtime=float(merged["max_runtime"]) if merged.get("max_runtime") is not None else None,
        target_count=int(merged["target_count"]) if merged.get("target_count") is not None else None,
        population=int(merged.get("population", 500)),
        iterations=int(merged.get("iterations", 1000)),
        mutation_rate=float(merged.get("mutation_rate", 0.05)),
        elite_fraction=float(merged.get("elite_fraction", 0.10)),
        experimental_signal=bool(getattr(args, "experimental_signal", False) or getattr(args, "experimental", False)),
        chemistry=str(getattr(args, "chemistry", "") or ""),
        plate_format=int(merged.get("plate_format", 96)),
    )
    assignment = AssignmentConfig(
        metric=_normalize_metric(str(merged.get("metric", "hamming"))),
        k=int(merged.get("k", 1)),
        reverse_complement_mode=str(merged.get("reverse_complement_mode", "warn")),
    )
    constraints = _constraints_from_args(args, length=options.length)
    constraints.alphabet = str(merged.get("alphabet", DNA))
    for key, attr in [
        ("min_hamming_distance", "min_hamming_distance"),
        ("min_levenshtein_distance", "min_levenshtein_distance"),
        ("min_sequence_levenshtein_distance", "min_sequence_levenshtein_distance"),
        ("gc_min", "gc_min"),
        ("gc_max", "gc_max"),
        ("max_homopolymer", "max_homopolymer"),
        ("max_dinucleotide_repeat", "max_dinucleotide_repeat"),
        ("avoid_reverse_complements", "avoid_reverse_complements"),
        ("avoid_self_complementarity", "avoid_self_complementarity"),
        ("avoid_prefix_collisions", "avoid_prefix_collisions"),
        ("avoid_suffix_collisions", "avoid_suffix_collisions"),
        ("avoid_ambiguous_bases", "avoid_ambiguous_bases"),
        ("cycle_balance_enabled", "cycle_balance_enabled"),
        ("cycle_min_base_fraction", "cycle_min_base_fraction"),
        ("cycle_max_base_fraction", "cycle_max_base_fraction"),
    ]:
        if key in merged:
            setattr(constraints, attr, merged[key])
    if "forbidden_motifs" in merged and isinstance(merged["forbidden_motifs"], list):
        constraints.forbidden_motifs = [str(item).upper() for item in merged["forbidden_motifs"]]
    return options, assignment, constraints


def _load_design_spec(path: Path) -> dict[str, object]:
    if not path.exists():
        raise ValueError(f"spec file not found: {path}")
    text = path.read_text(encoding="utf-8")
    if path.suffix.lower() == ".json":
        data = json.loads(text)
        return data if isinstance(data, dict) else {}
    return _parse_simple_yaml(text)


def _parse_simple_yaml(text: str) -> dict[str, object]:
    root: dict[str, object] = {}
    stack: list[tuple[int, dict[str, object] | list[object]]] = [(-1, root)]
    last_key_at_indent: dict[int, tuple[dict[str, object], str]] = {}
    for raw in text.splitlines():
        if not raw.strip() or raw.lstrip().startswith("#"):
            continue
        indent = len(raw) - len(raw.lstrip(" "))
        line = raw.strip()
        while stack and indent <= stack[-1][0]:
            stack.pop()
        parent = stack[-1][1]
        if line.startswith("- "):
            item = _parse_scalar(line[2:].strip())
            if not isinstance(parent, list):
                holder, key = last_key_at_indent[indent]
                new_list: list[object] = []
                holder[key] = new_list
                parent = new_list
                stack.append((indent, parent))
            assert isinstance(parent, list)
            parent.append(item)
            continue
        key, value = line.split(":", 1)
        key = key.strip()
        value = value.strip()
        if value == "":
            new_map: dict[str, object] = {}
            assert isinstance(parent, dict)
            parent[key] = new_map
            last_key_at_indent[indent + 2] = (parent, key)
            stack.append((indent, new_map))
        else:
            assert isinstance(parent, dict)
            parent[key] = _parse_scalar(value)
            last_key_at_indent[indent + 2] = (parent, key)
    return root


def _parse_scalar(value: str) -> object:
    value = value.strip()
    if (value.startswith('"') and value.endswith('"')) or (value.startswith("'") and value.endswith("'")):
        return value[1:-1]
    if value.lower() in {"true", "false"}:
        return value.lower() == "true"
    if value.startswith("[") and value.endswith("]"):
        inner = value[1:-1].strip()
        if not inner:
            return []
        return [_parse_scalar(part.strip()) for part in inner.split(",")]
    try:
        if any(ch in value for ch in [".", "e", "E"]):
            return float(value)
        return int(value)
    except ValueError:
        return value


def _as_dict(value: object) -> dict[str, object]:
    return value if isinstance(value, dict) else {}


def _normalize_metric(metric: str) -> str:
    if metric == "seqlev":
        return "sequence-levenshtein"
    return metric


def assign_query(query: str, records: Sequence[PanelRecord], metric: str, k: int) -> AssignmentOutcome:
    query = query.upper()
    if any(base not in DNA for base in query):
        return AssignmentOutcome(-1, -1, -1, 0, "invalid")
    metric = _normalize_metric(metric)
    matches: list[tuple[int, int]] = []
    distances: list[int] = []
    for i, record in enumerate(records):
        if metric == "hamming" and len(query) != len(record.sequence):
            continue
        dist = metric_distance(query, record.sequence, metric)
        distances.append(dist)
        if dist <= k:
            matches.append((i, dist))
    if not matches:
        return AssignmentOutcome(-1, -1, -1, 0, "none")
    best = min(dist for _i, dist in matches)
    best_matches = [(i, dist) for i, dist in matches if dist == best]
    second = min((dist for _i, dist in matches if dist > best), default=-1)
    if len(best_matches) == 1:
        return AssignmentOutcome(best_matches[0][0], best, second, len(matches), "unique")
    return AssignmentOutcome(-1, best, second, len(best_matches), "ambiguous")


def metric_distance(a: str, b: str, metric: str, constraints: PanelConstraints | None = None) -> int:
    metric = _normalize_metric(metric)
    if metric == "hamming":
        return hamming_distance(a, b)
    if metric == "sequence-levenshtein":
        constraints = constraints or PanelConstraints()
        left = constraints.left_flank[-constraints.context_window :] if constraints.context_window and constraints.left_flank else constraints.left_flank
        right = constraints.right_flank[: constraints.context_window] if constraints.context_window and constraints.right_flank else constraints.right_flank
        return levenshtein_distance(left + a + right, left + b + right)
    return levenshtein_distance(a, b)


def hamming_distance(a: str, b: str) -> int:
    if len(a) != len(b):
        return max(len(a), len(b))
    return sum(1 for x, y in zip(a, b) if x != y)


def levenshtein_distance(a: str, b: str) -> int:
    if a == b:
        return 0
    if len(a) < len(b):
        a, b = b, a
    previous = list(range(len(b) + 1))
    for i, ca in enumerate(a, start=1):
        current = [i]
        for j, cb in enumerate(b, start=1):
            current.append(min(previous[j] + 1, current[j - 1] + 1, previous[j - 1] + (ca != cb)))
        previous = current
    return previous[-1]


def variants_within(seq: str, metric: str, k: int) -> set[str]:
    metric = _normalize_metric(metric)
    variants = {seq}
    if k <= 0:
        return variants
    frontier = {seq}
    for _step in range(k):
        next_frontier: set[str] = set()
        for item in frontier:
            next_frontier.update(_one_edit_variants(item, include_indels=metric != "hamming"))
        variants.update(next_frontier)
        frontier = next_frontier
    return variants


def _one_edit_variants(seq: str, *, include_indels: bool) -> set[str]:
    variants: set[str] = set()
    for i, base in enumerate(seq):
        for alt in DNA:
            if alt != base:
                variants.add(seq[:i] + alt + seq[i + 1 :])
    if include_indels:
        for i in range(len(seq)):
            variants.add(seq[:i] + seq[i + 1 :])
        for i in range(len(seq) + 1):
            for base in DNA:
                variants.add(seq[:i] + base + seq[i:])
    return variants


def _assignment_risk_rows(records: Sequence[PanelRecord], assignment: AssignmentConfig, constraints: PanelConstraints) -> list[dict[str, object]]:
    if assignment.k < 0:
        return []
    if assignment.k > MAX_EXACT_CERTIFICATE_RADIUS:
        raise ValueError(f"exact safety certificate supports k <= {MAX_EXACT_CERTIFICATE_RADIUS}; refusing to certify k={assignment.k}")
    rows: list[dict[str, object]] = []
    enumerate_k = assignment.k
    for source_index, source in enumerate(records):
        for variant in sorted(variants_within(source.sequence, assignment.metric, enumerate_k)):
            outcome = assign_query(variant, records, assignment.metric, assignment.k)
            if outcome.status == "ambiguous":
                candidates = _best_candidate_ids(variant, records, assignment.metric, assignment.k)
                rows.append(
                    {
                        "source_barcode_id": source.barcode_id,
                        "source_sequence": source.sequence,
                        "variant": variant,
                        "risk_type": "ambiguous",
                        "assigned_barcode_id": "",
                        "assigned_sequence": "",
                        "best_distance": outcome.best_distance,
                        "match_count": outcome.match_count,
                        "candidate_ids": ",".join(candidates),
                        "metric": _normalize_metric(assignment.metric),
                        "k": assignment.k,
                        "certified_command": _certified_command(Path("barcodes.tsv"), assignment, constraints.length or len(source.sequence)),
                    }
                )
            elif outcome.status == "unique" and outcome.target_index != source_index:
                assigned = records[outcome.target_index]
                rows.append(
                    {
                        "source_barcode_id": source.barcode_id,
                        "source_sequence": source.sequence,
                        "variant": variant,
                        "risk_type": "silent_assignment",
                        "assigned_barcode_id": assigned.barcode_id,
                        "assigned_sequence": assigned.sequence,
                        "best_distance": outcome.best_distance,
                        "match_count": outcome.match_count,
                        "candidate_ids": assigned.barcode_id,
                        "metric": _normalize_metric(assignment.metric),
                        "k": assignment.k,
                        "certified_command": _certified_command(Path("barcodes.tsv"), assignment, constraints.length or len(source.sequence)),
                    }
                )
    return rows


def _best_candidate_ids(query: str, records: Sequence[PanelRecord], metric: str, k: int) -> list[str]:
    distances = [(record.barcode_id, metric_distance(query, record.sequence, metric)) for record in records if metric != "hamming" or len(query) == len(record.sequence)]
    within = [(barcode_id, dist) for barcode_id, dist in distances if dist <= k]
    if not within:
        return []
    best = min(dist for _barcode_id, dist in within)
    return [barcode_id for barcode_id, dist in within if dist == best]


def _safe_for_radius(records: Sequence[PanelRecord], metric: str, k: int, constraints: PanelConstraints) -> bool:
    assignment = AssignmentConfig(metric=metric, k=k)
    return not _assignment_risk_rows(records, assignment, constraints)


def _safe_by_distance_rule(records: Sequence[PanelRecord], metric: str, k: int, constraints: PanelConstraints) -> bool:
    min_dist = _minimum_pairwise_distance(records, metric, constraints)
    if min_dist is None:
        return True
    return min_dist >= 2 * k + 1


def pair_compatible(a: str, b: str, constraints: PanelConstraints) -> bool:
    if constraints.min_hamming_distance and hamming_distance(a, b) < constraints.min_hamming_distance:
        return False
    if constraints.min_levenshtein_distance and levenshtein_distance(a, b) < constraints.min_levenshtein_distance:
        return False
    if constraints.min_sequence_levenshtein_distance and metric_distance(a, b, "sequence-levenshtein", constraints) < constraints.min_sequence_levenshtein_distance:
        return False
    if constraints.avoid_reverse_complements:
        rc_b = reverse_complement(b)
        if a == rc_b:
            return False
        if len(a) == len(rc_b) and constraints.min_hamming_distance and hamming_distance(a, rc_b) < constraints.min_hamming_distance:
            return False
    if constraints.avoid_prefix_collisions and (a.startswith(b) or b.startswith(a)):
        return False
    if constraints.avoid_suffix_collisions and (a.endswith(b) or b.endswith(a)):
        return False
    return True


def candidate_filter_reason(seq: str, constraints: PanelConstraints) -> str:
    seq = seq.upper()
    if constraints.length is not None and len(seq) != constraints.length:
        return "wrong_length"
    if constraints.avoid_ambiguous_bases and any(base not in DNA for base in seq):
        return "ambiguous_base"
    if not constraints.allow_ambiguous_literals and any(base not in constraints.alphabet for base in seq):
        return "outside_alphabet"
    content = sequence_content(seq, constraints)
    if content["status"] == "fail":
        return ";".join(content["warnings"]) or "sequence_constraint"
    if constraints.left_flank or constraints.right_flank:
        context = context_risk(seq, constraints)
        if context:
            return "context_risk"
    return "pass"


def sequence_content(seq: str, constraints: PanelConstraints) -> dict[str, object]:
    gc = (seq.count("G") + seq.count("C")) / len(seq) if seq else 0.0
    homopolymer = max_homopolymer(seq)
    dinuc = max_dinucleotide_repeat(seq)
    motifs = [motif for motif in constraints.forbidden_motifs if motif and motif in seq]
    warnings: list[str] = []
    status = "pass"
    if constraints.avoid_ambiguous_bases and any(base not in DNA for base in seq):
        warnings.append("ambiguous_base")
        status = "fail"
    if gc < constraints.gc_min or gc > constraints.gc_max:
        warnings.append("gc_out_of_range")
        status = "fail"
    if homopolymer > constraints.max_homopolymer:
        warnings.append("homopolymer_too_long")
        status = "fail"
    if dinuc > constraints.max_dinucleotide_repeat:
        warnings.append("dinucleotide_repeat_too_long")
        status = "fail"
    if motifs:
        warnings.append("forbidden_motif")
        status = "fail"
    if constraints.avoid_self_complementarity and self_complement_score(seq) >= 0.75:
        warnings.append("self_complementarity")
        status = "fail"
    if len(set(seq)) <= 1:
        warnings.append("low_complexity")
        status = "fail"
    return {
        "gc": gc,
        "homopolymer_max": homopolymer,
        "dinucleotide_repeat_max": dinuc,
        "forbidden_motifs": ",".join(motifs),
        "self_complement_score": self_complement_score(seq),
        "warnings": warnings,
        "status": status,
    }


def max_homopolymer(seq: str) -> int:
    best = run = 0
    prev = ""
    for base in seq:
        run = run + 1 if base == prev else 1
        best = max(best, run)
        prev = base
    return best


def max_dinucleotide_repeat(seq: str) -> int:
    best = 0
    for i in range(max(0, len(seq) - 1)):
        motif = seq[i : i + 2]
        if len(motif) < 2:
            continue
        count = 0
        pos = i
        while seq[pos : pos + 2] == motif:
            count += 1
            pos += 2
        best = max(best, count)
    return best


def reverse_complement(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def self_complement_score(seq: str) -> float:
    if not seq:
        return 0.0
    rc = reverse_complement(seq)
    return sum(1 for a, b in zip(seq, rc) if a == b) / len(seq)


def context_risk(seq: str, constraints: PanelConstraints) -> list[str]:
    risks: list[str] = []
    left_join = constraints.left_flank + seq
    right_join = seq + constraints.right_flank
    for joined, name in [(left_join, "left_flank"), (right_join, "right_flank")]:
        if joined and max_homopolymer(joined) > constraints.max_homopolymer:
            risks.append(f"{name}_homopolymer")
        for motif in constraints.forbidden_motifs:
            if motif and motif in joined:
                risks.append(f"{name}_forbidden_motif:{motif}")
    return sorted(set(risks))


def individual_candidate_score(seq: str, constraints: PanelConstraints) -> float:
    content = sequence_content(seq, constraints)
    gc = float(content["gc"])
    gc_score = 1.0 - abs(gc - 0.5) * 2
    complexity = len(set(seq)) / 4
    homopolymer_score = max(0.0, 1.0 - (int(content["homopolymer_max"]) - 1) / max(1, constraints.max_homopolymer))
    return gc_score * 4 + complexity * 3 + homopolymer_score * 2 - len(content["warnings"]) * 5


def incremental_candidate_score(
    candidate: PanelRecord,
    selected: Sequence[PanelRecord],
    constraints: PanelConstraints,
    objective: str,
) -> float:
    weights = OBJECTIVE_WEIGHTS.get(objective, OBJECTIVE_WEIGHTS["strict-demux"])
    base = individual_candidate_score(candidate.sequence, constraints)
    if not selected:
        return base
    h = min(hamming_distance(candidate.sequence, item.sequence) for item in selected)
    lev = min(levenshtein_distance(candidate.sequence, item.sequence) for item in selected)
    seqlev = min(metric_distance(candidate.sequence, item.sequence, "sequence-levenshtein", constraints) for item in selected)
    rc_penalty = sum(1 for item in selected if candidate.sequence == reverse_complement(item.sequence))
    prefix_penalty = sum(1 for item in selected if candidate.sequence.startswith(item.sequence) or item.sequence.startswith(candidate.sequence))
    return (
        base
        + weights["min_pairwise_hamming_distance"] * h
        + weights["min_pairwise_levenshtein_distance"] * lev
        + weights["min_pairwise_sequence_levenshtein_distance"] * seqlev
        + weights["reverse_complement_collision_count"] * rc_penalty
        + weights["prefix_collision_count"] * prefix_penalty
    )


def panel_objective(records: Sequence[PanelRecord], constraints: PanelConstraints, objective: str) -> float:
    if not records:
        return 0.0
    weights = OBJECTIVE_WEIGHTS.get(objective, OBJECTIVE_WEIGHTS["strict-demux"])
    hamming_distances = _pairwise_distances(records, "hamming", constraints)
    lev_distances = _pairwise_distances(records, "levenshtein", constraints)
    seqlev_distances = _pairwise_distances(records, "sequence-levenshtein", constraints)
    min_h = min(hamming_distances) if hamming_distances else len(records[0].sequence)
    min_l = min(lev_distances) if lev_distances else len(records[0].sequence)
    min_s = min(seqlev_distances) if seqlev_distances else len(records[0].sequence)
    percentile = _percentile(hamming_distances, 5) if hamming_distances else min_h
    cycle_score = _cycle_balance_score(records, constraints)
    gc_values = [(record.sequence.count("G") + record.sequence.count("C")) / len(record.sequence) for record in records]
    gc_distribution = 1.0 - min(1.0, abs((sum(gc_values) / len(gc_values)) - 0.5) * 2)
    rc_collisions = sum(1 for a, b in combinations(records, 2) if a.sequence == reverse_complement(b.sequence))
    prefix_collisions = sum(1 for a, b in combinations(records, 2) if a.sequence.startswith(b.sequence) or b.sequence.startswith(a.sequence))
    forbidden_hits = sum(1 for record in records for motif in constraints.forbidden_motifs if motif and motif in record.sequence)
    low_complexity = sum(1 for record in records if len(set(record.sequence)) <= 2)
    score = 0.0
    score += weights["min_pairwise_hamming_distance"] * min_h
    score += weights["min_pairwise_levenshtein_distance"] * min_l
    score += weights["min_pairwise_sequence_levenshtein_distance"] * min_s
    score += weights["pairwise_distance_5th_percentile"] * percentile
    score += weights["edit_sphere_separation_score"] * max(0, min_h - 2)
    score += weights["cycle_balance_score"] * cycle_score
    score += weights["plate_layout_score"] * 1.0
    score += weights["gc_distribution_score"] * gc_distribution
    score += weights["context_separation_score"] * (1.0 if not any(context_risk(record.sequence, constraints) for record in records) else 0.0)
    score += weights["ambiguous_variant_count"] * 0
    score += weights["reverse_complement_collision_count"] * rc_collisions
    score += weights["prefix_collision_count"] * prefix_collisions
    score += weights["forbidden_motif_count"] * forbidden_hits
    score += weights["low_complexity_penalty"] * low_complexity
    return score


def _pairwise_rows(records: Sequence[PanelRecord], constraints: PanelConstraints) -> tuple[list[dict[str, object]], dict[str, dict[str, object]]]:
    rows: list[dict[str, object]] = []
    nearest: dict[str, dict[str, object]] = {record.barcode_id: {"distance": None, "neighbor": ""} for record in records}
    for a, b in combinations(records, 2):
        h = hamming_distance(a.sequence, b.sequence)
        lev = levenshtein_distance(a.sequence, b.sequence)
        seqlev = metric_distance(a.sequence, b.sequence, "sequence-levenshtein", constraints)
        rc_dist = hamming_distance(a.sequence, reverse_complement(b.sequence)) if len(a.sequence) == len(b.sequence) else levenshtein_distance(a.sequence, reverse_complement(b.sequence))
        warnings = []
        status = "pass"
        if h == 0 or lev == 0:
            warnings.append("duplicate")
            status = "fail"
        if h < constraints.min_hamming_distance:
            warnings.append("below_min_hamming")
            status = "warn" if status == "pass" else status
        if lev < constraints.min_levenshtein_distance:
            warnings.append("below_min_levenshtein")
            status = "warn" if status == "pass" else status
        if seqlev < constraints.min_sequence_levenshtein_distance:
            warnings.append("below_min_sequence_levenshtein")
            status = "warn" if status == "pass" else status
        if a.sequence == reverse_complement(b.sequence):
            warnings.append("reverse_complement_exact")
            status = "warn" if status == "pass" else status
        rows.append(
            {
                "barcode_id": a.barcode_id,
                "sequence": a.sequence,
                "other_id": b.barcode_id,
                "other_sequence": b.sequence,
                "hamming_distance": h,
                "levenshtein_distance": lev,
                "sequence_levenshtein_distance": seqlev,
                "reverse_complement_distance": rc_dist,
                "status": status,
                "warnings": ";".join(warnings),
            }
        )
        for record, other in [(a, b), (b, a)]:
            current = nearest[record.barcode_id]["distance"]
            if current is None or h < int(current):
                nearest[record.barcode_id] = {"distance": h, "neighbor": other.barcode_id}
    rows.sort(key=lambda row: (int(row["hamming_distance"]), int(row["levenshtein_distance"]), row["barcode_id"], row["other_id"]))
    return rows, nearest


def _target_safety_rows(
    records: Sequence[PanelRecord],
    content_rows: Sequence[dict[str, object]],
    risk_rows: Sequence[dict[str, object]],
) -> list[dict[str, object]]:
    content_by_id = {str(row["barcode_id"]): row for row in content_rows}
    risk_by_id: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in risk_rows:
        risk_by_id[str(row["source_barcode_id"])].append(row)
    rows: list[dict[str, object]] = []
    for record in records:
        risks = risk_by_id.get(record.barcode_id, [])
        content = content_by_id.get(record.barcode_id, {})
        status = "fail" if risks or content.get("status") == "fail" else "pass"
        rows.append(
            {
                "barcode_id": record.barcode_id,
                "sequence": record.sequence,
                "status": status,
                "ambiguous_variants": sum(1 for row in risks if row["risk_type"] == "ambiguous"),
                "silent_assignment_variants": sum(1 for row in risks if row["risk_type"] == "silent_assignment"),
                "nearest_hamming_neighbor": content.get("nearest_hamming_neighbor", ""),
                "gc": content.get("gc", ""),
                "homopolymer_max": content.get("homopolymer_max", ""),
                "warnings": content.get("warnings", ""),
            }
        )
    return rows


def _prefix_suffix_rows(records: Sequence[PanelRecord], assignment: AssignmentConfig, *, prefix: bool) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for a, b in combinations(records, 2):
        if prefix:
            collision = a.sequence.startswith(b.sequence) or b.sequence.startswith(a.sequence)
            overlap = _shared_prefix(a.sequence, b.sequence)
        else:
            collision = a.sequence.endswith(b.sequence) or b.sequence.endswith(a.sequence)
            overlap = _shared_suffix(a.sequence, b.sequence)
        near_overlap = overlap >= min(len(a.sequence), len(b.sequence)) - max(1, assignment.k)
        if collision or (len(a.sequence) != len(b.sequence) and near_overlap):
            rows.append(
                {
                    "barcode_id": a.barcode_id,
                    "sequence": a.sequence,
                    "other_id": b.barcode_id,
                    "other_sequence": b.sequence,
                    "overlap_length": overlap,
                    "risk": "exact_overlap" if collision else "near_overlap",
                    "metric": _normalize_metric(assignment.metric),
                    "k": assignment.k,
                }
            )
    return rows


def _reverse_complement_rows(records: Sequence[PanelRecord], assignment: AssignmentConfig, constraints: PanelConstraints) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for a, b in combinations(records, 2):
        rc_b = reverse_complement(b.sequence)
        dist = metric_distance(a.sequence, rc_b, assignment.metric, constraints)
        exact = a.sequence == rc_b
        if exact or dist <= max(assignment.k, 1):
            rows.append(
                {
                    "barcode_id": a.barcode_id,
                    "sequence": a.sequence,
                    "other_id": b.barcode_id,
                    "other_sequence": b.sequence,
                    "reverse_complement": rc_b,
                    "distance": dist,
                    "mode": assignment.reverse_complement_mode,
                    "risk": "exact_rc_collision" if exact else "near_rc_collision",
                }
            )
    return rows


def _content_rows(records: Sequence[PanelRecord], constraints: PanelConstraints, nearest: dict[str, dict[str, object]]) -> list[dict[str, object]]:
    rows = []
    for record in records:
        content = sequence_content(record.sequence, constraints)
        rows.append(
            {
                "barcode_id": record.barcode_id,
                "sequence": record.sequence,
                "length": len(record.sequence),
                "gc": f"{float(content['gc']):.8f}",
                "homopolymer_max": content["homopolymer_max"],
                "dinucleotide_repeat_max": content["dinucleotide_repeat_max"],
                "self_complement_score": f"{float(content['self_complement_score']):.8f}",
                "forbidden_motifs": content["forbidden_motifs"],
                "nearest_hamming_neighbor": nearest.get(record.barcode_id, {}).get("neighbor", ""),
                "nearest_hamming_distance": nearest.get(record.barcode_id, {}).get("distance", ""),
                "status": content["status"],
                "warnings": ";".join(content["warnings"]),
            }
        )
    return rows


def _cycle_balance_rows(records: Sequence[PanelRecord], constraints: PanelConstraints) -> tuple[list[dict[str, object]], int, float]:
    if not records:
        return [], 0, 0.0
    max_len = max(len(record.sequence) for record in records)
    rows: list[dict[str, object]] = []
    warnings = 0
    score_total = 0.0
    for pos in range(max_len):
        counts = Counter(record.sequence[pos] for record in records if pos < len(record.sequence))
        total = sum(counts.values())
        fractions = {base: counts.get(base, 0) / total if total else 0.0 for base in DNA}
        status = "pass"
        if constraints.cycle_balance_enabled and total >= 8:
            if any(frac < constraints.cycle_min_base_fraction or frac > constraints.cycle_max_base_fraction for frac in fractions.values()):
                status = "warn"
                warnings += 1
        elif constraints.cycle_balance_enabled and total < 8:
            status = "warn"
            warnings += 1
        spread = max(fractions.values()) - min(fractions.values()) if fractions else 1.0
        score_total += max(0.0, 1.0 - spread)
        rows.append(
            {
                "cycle": pos + 1,
                "a_fraction": f"{fractions['A']:.8f}",
                "c_fraction": f"{fractions['C']:.8f}",
                "g_fraction": f"{fractions['G']:.8f}",
                "t_fraction": f"{fractions['T']:.8f}",
                "status": status,
            }
        )
    return rows, warnings, score_total / len(rows) if rows else 0.0


def _context_rows(records: Sequence[PanelRecord], constraints: PanelConstraints) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for record in records:
        risks = context_risk(record.sequence, constraints)
        if risks or constraints.left_flank or constraints.right_flank:
            rows.append(
                {
                    "barcode_id": record.barcode_id,
                    "sequence": record.sequence,
                    "left_flank": constraints.left_flank,
                    "right_flank": constraints.right_flank,
                    "flanked_sequence": constraints.left_flank + record.sequence + constraints.right_flank,
                    "context_window": constraints.context_window,
                    "status": "warn" if risks else "pass",
                    "risks": ";".join(risks),
                }
            )
    return rows


def _flanked_rows(records: Sequence[PanelRecord], constraints: PanelConstraints) -> list[dict[str, object]]:
    return [
        {
            "barcode_id": record.barcode_id,
            "sequence": record.sequence,
            "left_flank": constraints.left_flank,
            "right_flank": constraints.right_flank,
            "flanked_sequence": constraints.left_flank + record.sequence + constraints.right_flank,
            "context_window": constraints.context_window,
            "cross_boundary_homopolymer": max_homopolymer(constraints.left_flank + record.sequence + constraints.right_flank),
            "risks": ";".join(context_risk(record.sequence, constraints)),
        }
        for record in records
    ]


def _minimum_pairwise_distance(records: Sequence[PanelRecord], metric: str, constraints: PanelConstraints | None = None) -> int | None:
    distances = _pairwise_distances(records, metric, constraints)
    return min(distances) if distances else None


def _nearest_distance(record: PanelRecord, others: Sequence[PanelRecord], metric: str, constraints: PanelConstraints | None = None) -> int | str:
    if not others:
        return ""
    return min(metric_distance(record.sequence, other.sequence, metric, constraints) for other in others)


def _pairwise_distances(records: Sequence[PanelRecord], metric: str, constraints: PanelConstraints | None = None) -> list[int]:
    return [metric_distance(a.sequence, b.sequence, metric, constraints) for a, b in combinations(records, 2)]


def _guaranteed_correction(min_distance: int | None) -> int:
    return max(0, (min_distance - 1) // 2) if min_distance is not None else 0


def _guaranteed_detection(min_distance: int | None) -> int:
    return max(0, min_distance - 1) if min_distance is not None else 0


def _panel_grade(
    status: str,
    records: Sequence[PanelRecord],
    assignment: AssignmentConfig,
    constraints: PanelConstraints,
    min_hamming: int | None,
    min_lev: int | None,
    ambiguous_rows: Sequence[dict[str, object]],
    silent_rows: Sequence[dict[str, object]],
    rc_rows: Sequence[dict[str, object]],
    prefix_rows: Sequence[dict[str, object]],
    suffix_rows: Sequence[dict[str, object]],
) -> str:
    if status == "fail":
        if ambiguous_rows or silent_rows or prefix_rows or suffix_rows:
            return "F"
        return "D"
    if not _safe_by_distance_rule(records, assignment.metric, 1, constraints):
        return "C"
    if status == "warn" or rc_rows:
        return "B"
    if min_hamming is not None and min_lev is not None and min_hamming >= constraints.min_hamming_distance + 1 and min_lev >= constraints.min_levenshtein_distance + 1:
        return "A+"
    return "A"


def _narrow_margin(records: Sequence[PanelRecord], assignment: AssignmentConfig) -> bool:
    min_dist = _minimum_pairwise_distance(records, assignment.metric)
    return bool(min_dist is not None and min_dist <= 2 * assignment.k + 1)


def _cycle_balance_score(records: Sequence[PanelRecord], constraints: PanelConstraints) -> float:
    _rows, _warnings, score = _cycle_balance_rows(records, constraints)
    return score


def _compatibility_degree(record: PanelRecord, candidates: Sequence[PanelRecord], constraints: PanelConstraints) -> int:
    return sum(1 for other in candidates if other != record and pair_compatible(record.sequence, other.sequence, constraints))


def _renumber_records(records: Sequence[PanelRecord], prefix: str = "BC") -> list[PanelRecord]:
    return [PanelRecord(f"{prefix}{i + 1:03d}", record.sequence) for i, record in enumerate(records)]


def mutate_sequence(seq: str, rng: random.Random, substitution_rate: float, insertion_rate: float, deletion_rate: float) -> str:
    out: list[str] = []
    for base in seq:
        if rng.random() < deletion_rate:
            continue
        if rng.random() < substitution_rate:
            choices = [item for item in DNA if item != base]
            out.append(rng.choice(choices))
        else:
            out.append(base)
        if rng.random() < insertion_rate:
            out.append(rng.choice(DNA))
    if rng.random() < insertion_rate:
        out.append(rng.choice(DNA))
    return "".join(out)


def _plate_shape(plate: int) -> tuple[int, int]:
    if plate == 24:
        return 4, 6
    if plate == 384:
        return 16, 24
    return 8, 12


def _plate_order(records: Sequence[PanelRecord]) -> list[PanelRecord]:
    if len(records) <= 2:
        return list(records)
    remaining = list(records)
    ordered = [remaining.pop(0)]
    while remaining:
        remaining.sort(key=lambda record: -min(hamming_distance(record.sequence, placed.sequence) for placed in ordered))
        ordered.append(remaining.pop(0))
    return ordered


def _neighbor_rows(layout_rows: Sequence[dict[str, object]], records: Sequence[PanelRecord]) -> list[dict[str, object]]:
    by_well = {str(row["well"]): row for row in layout_rows}
    rows: list[dict[str, object]] = []
    for row in layout_rows:
        well = str(row["well"])
        row_letter = str(row["row"])
        col = int(row["column"])
        for neighbor in [f"{row_letter}{col + 1}", f"{chr(ord(row_letter) + 1)}{col}" if len(row_letter) == 1 and row_letter < "Z" else ""]:
            if neighbor in by_well:
                other = by_well[neighbor]
                rows.append(
                    {
                        "well": well,
                        "neighbor_well": neighbor,
                        "barcode_id": row["barcode_id"],
                        "neighbor_barcode_id": other["barcode_id"],
                        "hamming_distance": hamming_distance(str(row["sequence"]), str(other["sequence"])),
                    }
                )
    return rows


def _write_picklist(path: Path, layout_rows: Sequence[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["Source Well", "Sample ID", "Barcode Sequence"])
        writer.writeheader()
        for row in layout_rows:
            writer.writerow({"Source Well": row["well"], "Sample ID": row["barcode_id"], "Barcode Sequence": row["sequence"]})


def _write_plate_svg(path: Path, layout_rows: Sequence[dict[str, object]], rows: int, cols: int) -> None:
    cell = 36
    margin = 28
    by_well = {str(row["well"]): row for row in layout_rows}
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{cols * cell + margin * 2}" height="{rows * cell + margin * 2}" viewBox="0 0 {cols * cell + margin * 2} {rows * cell + margin * 2}">',
        '<rect width="100%" height="100%" fill="white"/>',
    ]
    for r in range(rows):
        for c in range(cols):
            well = f"{chr(ord('A') + r)}{c + 1}"
            x = margin + c * cell
            y = margin + r * cell
            fill = "#dff3ea" if well in by_well else "#f4f6f5"
            parts.append(f'<rect x="{x}" y="{y}" width="30" height="30" rx="4" fill="{fill}" stroke="#557066"/>')
            parts.append(f'<text x="{x + 15}" y="{y + 19}" font-size="9" text-anchor="middle" fill="#10201a">{html.escape(well)}</text>')
    parts.append("</svg>")
    path.write_text("\n".join(parts) + "\n", encoding="utf-8")


def _write_illumina_samplesheet(path: Path, records: Sequence[PanelRecord]) -> None:
    lines = ["[Header]", "IEMFileVersion,5", "Workflow,GenerateFASTQ", "", "[Reads]", "151", "151", "", "[Data]", "Sample_ID,index"]
    for record in records:
        lines.append(f"{record.barcode_id},{record.sequence}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_lab_readme(path: Path, records: Sequence[PanelRecord], assignment: AssignmentConfig, source_panel: Path) -> None:
    length = _single_length(records) or "auto"
    command = _certified_command(source_panel, assignment, length)
    text = "\n".join(
        [
            "# DotMatch Barcode Panel Lab Notes",
            "",
            "Every DotMatch barcode panel ships with a machine-checkable safety certificate.",
            "",
            "Use this panel for known-target assignment under the certified settings below.",
            "",
            "```bash",
            command,
            "```",
            "",
            "Scope: barcode design for known-target assignment, not general genome alignment, not UMI entropy generation, not basecalling.",
            "",
        ]
    )
    path.write_text(text, encoding="utf-8")


def _write_signal_outputs(out_dir: Path, records: Sequence[PanelRecord], options: DesignOptions, certificate: dict[str, object]) -> None:
    templates = []
    for record in records:
        signal = "-".join(str(_base_signal_value(base, options.chemistry)) for base in record.sequence)
        templates.append({"barcode_id": record.barcode_id, "sequence": record.sequence, "predicted_template": signal})
    _write_tsv(out_dir / "predicted_signal_templates.tsv", ["barcode_id", "sequence", "predicted_template"], templates)
    matrix_rows = []
    confusable = []
    for a, b in combinations(records, 2):
        dist = sum(abs(_base_signal_value(x, options.chemistry) - _base_signal_value(y, options.chemistry)) for x, y in zip(a.sequence, b.sequence))
        matrix_rows.append({"barcode_id": a.barcode_id, "other_id": b.barcode_id, "signal_distance": f"{dist:.4f}"})
        confusable.append((dist, a.barcode_id, b.barcode_id))
    _write_tsv(out_dir / "signal_distance_matrix.tsv", ["barcode_id", "other_id", "signal_distance"], matrix_rows)
    confusable_rows = [
        {"barcode_id": a, "other_id": b, "signal_distance": f"{dist:.4f}"}
        for dist, a, b in sorted(confusable)[:20]
    ]
    _write_tsv(out_dir / "top_confusable_signal_pairs.tsv", ["barcode_id", "other_id", "signal_distance"], confusable_rows)
    _write_json(
        out_dir / "signal_design_report.json",
        {
            "schema_version": PANEL_SCHEMA_VERSION,
            "experimental": True,
            "chemistry": options.chemistry or "unspecified",
            "symbolic_safety_status": certificate["status"],
            "hard_rule": "Signal design cannot override symbolic safety.",
            "median_signal_margin": _median_float([dist for dist, _a, _b in confusable]),
            "worst_pair_signal_margin": min((dist for dist, _a, _b in confusable), default=None),
            "top_1_signal_accuracy": None,
            "top_3_recall": None,
            "top_5_recall": None,
            "false_exclusion_rate": None,
            "latency_per_read": None,
        },
    )


def _base_signal_value(base: str, chemistry: str) -> float:
    base_values = {"A": 1.0, "C": 2.1, "G": 3.0, "T": 4.2}
    offset = 0.2 if chemistry else 0.0
    return base_values.get(base, 0.0) + offset


def _write_panel_report(
    path: Path,
    summary: dict[str, object],
    records: Sequence[PanelRecord],
    pair_rows: Sequence[dict[str, object]],
    ambiguous_rows: Sequence[dict[str, object]],
    rc_rows: Sequence[dict[str, object]],
    prefix_rows: Sequence[dict[str, object]],
    suffix_rows: Sequence[dict[str, object]],
) -> None:
    closest = list(pair_rows)[:10]
    closest_rows = "".join(
        f"<tr><td>{html.escape(str(row['barcode_id']))}</td><td>{html.escape(str(row['other_id']))}</td><td>{html.escape(str(row['hamming_distance']))}</td><td>{html.escape(str(row['levenshtein_distance']))}</td><td>{html.escape(str(row['warnings']))}</td></tr>"
        for row in closest
    )
    warning = ""
    if summary["status"] == "fail":
        warning = (
            f"<div class=\"warning\"><strong>Do not use --metric {html.escape(str(summary['assignment_metric']))} --k {html.escape(str(summary['configured_assignment_k']))} with this panel.</strong> "
            f"Reason: {html.escape(str(summary['ambiguous_error_spheres']))} ambiguous variants and {html.escape(str(summary['silent_assignment_risk']))} silent assignment risks.</div>"
        )
    command = html.escape(str(summary["certified_dotmatch_command"]))
    document = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <title>DotMatch Panel Safety Report</title>
  <style>
    body {{ margin: 2rem; color: #111816; font-family: system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; }}
    table {{ border-collapse: collapse; width: 100%; }}
    th, td {{ border: 1px solid #d9e3df; padding: 0.45rem; text-align: left; }}
    th {{ background: #eef6f2; }}
    code, pre {{ background: #eef6f2; padding: 0.25rem; }}
    .warning {{ border-left: 4px solid #b3261e; background: #fff2f1; padding: 0.75rem; margin: 1rem 0; }}
  </style>
</head>
<body>
  <h1>DotMatch Panel Safety Report</h1>
  <p>Status: <strong>{html.escape(str(summary['status']).upper())}</strong>. Panel grade: <strong>{html.escape(str(summary['panel_grade']))}</strong>.</p>
  {warning}
  <h2>1. Can I use this panel?</h2>
  <p>{html.escape('Yes' if summary['status'] != 'fail' else 'No under the configured correction radius.')}</p>
  <h2>2. Safe settings</h2>
  <ul>
    <li>k=0 exact: {html.escape(str(summary['safe_for_k0']))}</li>
    <li>k=1 hamming: {html.escape(str(summary['safe_for_k1_hamming']))}</li>
    <li>k=1 levenshtein: {html.escape(str(summary['safe_for_k1_levenshtein']))}</li>
    <li>k=2 hamming: {html.escape(str(summary['safe_for_k2_hamming']))}</li>
    <li>k=2 levenshtein: {html.escape(str(summary['safe_for_k2_levenshtein']))}</li>
  </ul>
  <h2>3. Unsafe settings</h2>
  <p>Configured metric={html.escape(str(summary['assignment_metric']))}, k={html.escape(str(summary['configured_assignment_k']))}, ambiguous variants={html.escape(str(summary['ambiguous_error_spheres']))}, silent assignment risks={html.escape(str(summary['silent_assignment_risk']))}.</p>
  <h2>4. Closest barcode pairs</h2>
  <table><thead><tr><th>Barcode</th><th>Other</th><th>Hamming</th><th>Levenshtein</th><th>Warnings</th></tr></thead><tbody>{closest_rows}</tbody></table>
  <h2>5. DotMatch command</h2>
  <pre>{command}</pre>
  <h2>Safety Certificate</h2>
  <p>Every DotMatch barcode panel ships with a machine-checkable safety certificate. This report is for known-target assignment, not general genome alignment, not UMI entropy generation, not basecalling.</p>
  <p>Reverse-complement warnings: {len(rc_rows)}. Prefix collisions: {len(prefix_rows)}. Suffix collisions: {len(suffix_rows)}. Barcodes: {len(records)}.</p>
</body>
</html>
"""
    path.write_text(document, encoding="utf-8")


def _write_basic_html(path: Path, title: str, lines: Sequence[str]) -> None:
    items = "".join(f"<li>{html.escape(line)}</li>" for line in lines)
    path.write_text(
        f"<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\" /><title>{html.escape(title)}</title></head><body><h1>{html.escape(title)}</h1><ul>{items}</ul></body></html>\n",
        encoding="utf-8",
    )


def _write_json(path: Path, data: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _write_tsv(path: Path, columns: Sequence[str], rows: Sequence[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(columns), delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def _copy_if_exists(source: Path, dest: Path) -> None:
    if source.exists():
        dest.write_text(source.read_text(encoding="utf-8"), encoding="utf-8")


def _stamp_command(rows: Sequence[dict[str, object]], command: str) -> None:
    for row in rows:
        row["certified_command"] = command


def _certified_command(panel_path: Path, assignment: AssignmentConfig, length: int | str) -> str:
    return (
        "dotmatch demux "
        f"--barcodes {panel_path} "
        "--reads pooled.fastq.gz "
        "--barcode-start 0 "
        f"--barcode-length {length} "
        f"--k {assignment.k} "
        f"--metric {_normalize_metric(assignment.metric)} "
        "--out-dir demuxed/"
    )


def _single_length(records: Sequence[PanelRecord]) -> int | None:
    lengths = {len(record.sequence) for record in records}
    return next(iter(lengths)) if len(lengths) == 1 else None


def _first_present(names: Sequence[str], candidates: Sequence[str]) -> str | None:
    by_lower = {name.lower(): name for name in names}
    for candidate in candidates:
        if candidate.lower() in by_lower:
            return by_lower[candidate.lower()]
    return None


def _shared_prefix(a: str, b: str) -> int:
    count = 0
    for x, y in zip(a, b):
        if x != y:
            break
        count += 1
    return count


def _shared_suffix(a: str, b: str) -> int:
    return _shared_prefix(a[::-1], b[::-1])


def _percentile(values: Sequence[int], percentile: float) -> float:
    if not values:
        return 0.0
    sorted_values = sorted(values)
    index = min(len(sorted_values) - 1, max(0, int(math.ceil((percentile / 100) * len(sorted_values))) - 1))
    return float(sorted_values[index])


def _median_float(values: Sequence[float]) -> float | None:
    if not values:
        return None
    sorted_values = sorted(values)
    mid = len(sorted_values) // 2
    if len(sorted_values) % 2:
        return sorted_values[mid]
    return (sorted_values[mid - 1] + sorted_values[mid]) / 2


def _target_safety_columns() -> list[str]:
    return [
        "barcode_id",
        "sequence",
        "status",
        "ambiguous_variants",
        "silent_assignment_variants",
        "nearest_hamming_neighbor",
        "gc",
        "homopolymer_max",
        "warnings",
        "certified_command",
    ]


def _pairwise_columns() -> list[str]:
    return [
        "barcode_id",
        "sequence",
        "other_id",
        "other_sequence",
        "hamming_distance",
        "levenshtein_distance",
        "sequence_levenshtein_distance",
        "reverse_complement_distance",
        "status",
        "warnings",
        "certified_command",
    ]


def _risk_columns() -> list[str]:
    return [
        "source_barcode_id",
        "source_sequence",
        "variant",
        "risk_type",
        "assigned_barcode_id",
        "assigned_sequence",
        "best_distance",
        "match_count",
        "candidate_ids",
        "metric",
        "k",
        "certified_command",
    ]


def _prefix_columns() -> list[str]:
    return ["barcode_id", "sequence", "other_id", "other_sequence", "overlap_length", "risk", "metric", "k", "certified_command"]


def _rc_columns() -> list[str]:
    return ["barcode_id", "sequence", "other_id", "other_sequence", "reverse_complement", "distance", "mode", "risk", "certified_command"]


def _content_columns() -> list[str]:
    return [
        "barcode_id",
        "sequence",
        "length",
        "gc",
        "homopolymer_max",
        "dinucleotide_repeat_max",
        "self_complement_score",
        "forbidden_motifs",
        "nearest_hamming_neighbor",
        "nearest_hamming_distance",
        "status",
        "warnings",
        "certified_command",
    ]


def _cycle_columns() -> list[str]:
    return ["cycle", "a_fraction", "c_fraction", "g_fraction", "t_fraction", "status", "certified_command"]


def _context_columns() -> list[str]:
    return ["barcode_id", "sequence", "left_flank", "right_flank", "flanked_sequence", "context_window", "status", "risks", "certified_command"]


def _flanked_columns() -> list[str]:
    return [
        "barcode_id",
        "sequence",
        "left_flank",
        "right_flank",
        "flanked_sequence",
        "context_window",
        "cross_boundary_homopolymer",
        "risks",
        "certified_command",
    ]


def _simulation_columns() -> list[str]:
    return [
        "read_id",
        "true_id",
        "true_sequence",
        "observed_sequence",
        "assigned_id",
        "assigned_sequence",
        "distance",
        "status",
        "false_assignment",
        "certified_command",
    ]


def _confusion_columns() -> list[str]:
    return ["barcode_id", "true_reads", "unique_correct", "recall", "precision", "most_common_outcome", "certified_command"]
