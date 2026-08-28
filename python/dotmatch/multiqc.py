"""
DotMatch MultiQC module / plugin helpers.

This provides a "real" MultiQC module (beyond the custom-content TSV example)
that can be used by dropping this file into MultiQC's module search path or
by contributing similar code to the main MultiQC repo.

Usage (advanced):
- PYTHONPATH=... multiqc ... --module dotmatch  (if registered)
- Or copy the DotMatchModule class into your own multiqc/modules/dotmatch/

It parses:
- sample_qc.tsv (and *sample_qc.tsv)
- summary.json (and *summary.json)
- crispr_qc.summary.tsv
- assay_manifest.summary.tsv
- panel_summary.json
- top_unmatched.tsv-style diagnostics
and adds summary tables and plots to the MultiQC report.

See examples/workflows/multiqc/ for the simple custom-content alternative
that works today with no extra code.
"""

from __future__ import annotations

import csv
import json
import logging
from pathlib import Path
from typing import Any

log = logging.getLogger(__name__)

# Lazy imports: only fail when actually using the MultiQC integration
# (users who want it do: pip install multiqc "dotmatch[multiqc]" or similar)
# This keeps the core package lightweight while allowing MultiQC to load the
# integration when the optional dependency is installed.
#
# Registration happens via pyproject.toml entry point:
#   [project.entry-points."multiqc.modules.v1"]
#   dotmatch = "dotmatch.multiqc:DotMatchModule"
# After `pip install "dotmatch[multiqc]"`, `multiqc ...` should auto-discover it
# (or use --module dotmatch explicitly).
_HAS_MULTIQC = False
try:
    from multiqc import BaseMultiqcModule  # type: ignore
    from multiqc import config as _mqc_config  # type: ignore
    from multiqc.plots import table as _mqc_table  # type: ignore
    _HAS_MULTIQC = True
except Exception:  # noqa: BLE001
    BaseMultiqcModule = object  # type: ignore
    _mqc_table = None  # type: ignore
    _mqc_config = None  # type: ignore


DOTMATCH_SEARCH_PATTERNS: dict[str, dict[str, Any]] = {
    "dotmatch/sample_qc": {"fn": "*sample_qc.tsv"},
    # Keep the specific panel pattern before the generic summary pattern.
    # MultiQC stops after the first non-shared filename match, including an
    # excluded generic match.
    "dotmatch/panel_summary": {"fn": "*panel_summary.json"},
    "dotmatch/summary_json": {
        "fn": "*summary.json",
        "exclude_fn": "*panel_summary.json",
    },
    "dotmatch/crispr_qc": {"fn": "*crispr_qc.summary.tsv"},
    "dotmatch/assay_manifest": {"fn": "*assay_manifest.summary.tsv"},
    "dotmatch/top_unmatched": {"fn": "*top_unmatched.tsv"},
}


INTEGER_FIELDS = {
    "ambiguous_error_spheres",
    "ambiguous_reads",
    "assigned_corrected",
    "assigned_exact",
    "assigned_reads",
    "assigned_unique",
    "candidates_verified",
    "collision_pairs",
    "configured_assignment_k",
    "exact_reads",
    "invalid",
    "invalid_reads",
    "k1_del_reads",
    "k1_ins_reads",
    "k1_rescued_reads",
    "k1_sub_reads",
    "minimum_hamming_distance",
    "n_barcodes",
    "no_match_reads",
    "sample_count",
    "silent_assignment_risk",
    "targets_observed",
    "total_count",
    "total_reads",
    "unmatched",
    "valid_extracted_reads",
    "warning_count",
    "zero_count_targets",
}

FLOAT_FIELDS = {
    "ambiguous_rate",
    "assignment_rate",
    "coverage_fraction",
    "exact_rate",
    "gini_index",
    "invalid_rate",
    "no_match_rate",
    "rescue_rate",
    "top_1pct_fraction",
    "top_1pct_read_fraction",
    "zero_count_fraction",
}


def load_config() -> None:
    """Register DotMatch search patterns before MultiQC indexes input files."""
    if _mqc_config is None:
        return
    search_patterns = getattr(_mqc_config, "sp", None)
    if not isinstance(search_patterns, dict):
        return
    for key, pattern in DOTMATCH_SEARCH_PATTERNS.items():
        search_patterns.setdefault(key, pattern)


def _get_table_plot():
    if not _HAS_MULTIQC or _mqc_table is None:
        raise ImportError("multiqc is required to use the DotMatch MultiQC module")
    return _mqc_table.plot


def _coerce_value(name: str, value: Any) -> Any:
    if value is None or value == "":
        return ""
    if name in INTEGER_FIELDS:
        try:
            return int(value)
        except (TypeError, ValueError):
            return value
    if name in FLOAT_FIELDS:
        try:
            return float(value)
        except (TypeError, ValueError):
            return value
    return value


def _read_tsv_rows(path: str | Path) -> list[dict[str, Any]]:
    path = Path(path)
    with path.open(encoding="utf-8", newline="") as fh:
        return [
            {key: _coerce_value(key, value) for key, value in row.items()}
            for row in csv.DictReader(fh, delimiter="\t")
            if any(value not in (None, "") for value in row.values())
        ]


def _multiqc_file_path(record: Any) -> Path:
    """Return a filesystem path from a MultiQC find_log_files record."""
    if isinstance(record, dict):
        filename = record.get("path") or record.get("fn")
        root = record.get("root") or record.get("dir") or ""
        if filename is None:
            raise ValueError(f"MultiQC file record has no filename: {record!r}")
        path = Path(str(filename))
        if path.is_absolute() or not root:
            return path
        return Path(root) / path
    return Path(record)


# Pure, dependency-free parsers. These can be used from Python code (notebooks,
# custom pipelines, or other report generators) while preserving the documented
# DotMatch semantics when consuming TSV artifacts outside MultiQC.
def parse_sample_qc_tsv(path: str | Path) -> dict[str, dict[str, Any]]:
    """Parse DotMatch sample_qc.tsv into dict of per-sample metrics.

    Columns and semantics exactly as documented in docs/schemas.md.
    assignment_rate etc. denominator is valid_extracted_reads; ambiguous reads
    are never counted toward any target (core scientific accuracy property).
    """
    path = Path(path)
    data: dict[str, dict[str, Any]] = {}
    for row in _read_tsv_rows(path):
        sname = str(row.get("sample_id") or row.get("sample") or path.stem)
        data[sname] = {
            "total_reads": row.get("total_reads", ""),
            "valid_extracted_reads": row.get("valid_extracted_reads", ""),
            "assigned_reads": row.get("assigned_reads", ""),
            "exact_reads": row.get("exact_reads", ""),
            "k1_rescued_reads": row.get("k1_rescued_reads", ""),
            "ambiguous_reads": row.get("ambiguous_reads", ""),
            "no_match_reads": row.get("no_match_reads", ""),
            "invalid_reads": row.get("invalid_reads", ""),
            "assignment_rate": row.get("assignment_rate", ""),
            "exact_rate": row.get("exact_rate", ""),
            "rescue_rate": row.get("rescue_rate", ""),
            "ambiguous_rate": row.get("ambiguous_rate", ""),
            "no_match_rate": row.get("no_match_rate", ""),
            "targets_observed": row.get("targets_observed", ""),
            "zero_count_targets": row.get("zero_count_targets", ""),
            "gini_index": row.get("gini_index", ""),
            "top_1pct_read_fraction": row.get("top_1pct_read_fraction", ""),
            "candidates_verified": row.get("candidates_verified", ""),
        }
    return data

def parse_crispr_qc_summary_tsv(path: str | Path) -> dict[str, dict[str, Any]]:
    """Parse crispr_qc.summary.tsv (guide library representation + QC)."""
    path = Path(path)
    data: dict[str, dict[str, Any]] = {}
    for row in _read_tsv_rows(path):
        sname = str(row.get("sample_id") or path.stem)
        data[sname] = {
            "qc_status": row.get("qc_status", ""),
            "total_count": row.get("total_count", ""),
            "coverage_fraction": row.get("coverage_fraction", ""),
            "zero_count_fraction": row.get("zero_count_fraction", ""),
            "gini_index": row.get("gini_index", ""),
            "top_1pct_fraction": row.get("top_1pct_fraction", ""),
            "assignment_rate": row.get("assignment_rate", ""),
            "ambiguous_rate": row.get("ambiguous_rate", ""),
            "no_match_rate": row.get("no_match_rate", ""),
            "invalid_rate": row.get("invalid_rate", ""),
        }
    return data

def parse_assay_manifest_summary_tsv(path: str | Path) -> dict[str, dict[str, Any]]:
    """Parse assay_manifest.summary.tsv for run provenance and links."""
    path = Path(path)
    data: dict[str, dict[str, Any]] = {}
    for row in _read_tsv_rows(path):
        key = f"{row.get('assay_type', 'assay')}_{row.get('mode', 'run')}"
        data[key] = {
            "mode": row.get("mode", ""),
            "assay_type": row.get("assay_type", ""),
            "status": row.get("status", ""),
            "native_version": row.get("native_version", ""),
            "sample_count": row.get("sample_count", ""),
            "autopsy_triggered": row.get("autopsy_triggered", ""),
            "warning_count": row.get("warning_count", ""),
            "primary_report": row.get("primary_report", ""),
        }
    return data


def parse_summary_json(path: str | Path) -> dict[str, dict[str, Any]]:
    """Parse DotMatch summary.json into a single MultiQC-friendly row."""
    path = Path(path)
    summary = json.loads(path.read_text(encoding="utf-8"))
    sample = summary.get("sample_label") or summary.get("sample") or path.stem.replace(".summary", "")
    samples = summary.get("samples")
    if isinstance(samples, list) and samples:
        sample = ",".join(str(item.get("sample_id", item.get("id", ""))) for item in samples if isinstance(item, dict)) or sample
    return {
        str(sample): {
            "total_reads": _coerce_value("total_reads", summary.get("total_reads", "")),
            "assigned_unique": _coerce_value("assigned_unique", summary.get("assigned_unique", summary.get("assigned_reads", ""))),
            "assigned_exact": _coerce_value("assigned_exact", summary.get("assigned_exact", summary.get("exact_reads", ""))),
            "assigned_corrected": _coerce_value("assigned_corrected", summary.get("assigned_corrected", "")),
            "ambiguous": _coerce_value("ambiguous", summary.get("ambiguous", summary.get("ambiguous_reads", ""))),
            "unmatched": _coerce_value("unmatched", summary.get("unmatched", summary.get("no_match_reads", ""))),
            "invalid": _coerce_value("invalid", summary.get("invalid", summary.get("invalid_reads", ""))),
            "assignment_rate": _coerce_value("assignment_rate", summary.get("assignment_rate", "")),
            "ambiguous_rate": _coerce_value("ambiguous_rate", summary.get("ambiguous_rate", "")),
            "assignment_engine": summary.get("assignment_engine", ""),
        }
    }


def parse_panel_summary_json(path: str | Path) -> dict[str, dict[str, Any]]:
    """Parse dotmatch panel check panel_summary.json into a MultiQC row."""
    path = Path(path)
    summary = json.loads(path.read_text(encoding="utf-8"))
    key = str(summary.get("panel_id") or path.parent.name or path.stem)
    return {
        key: {
            "status": summary.get("status", ""),
            "panel_grade": summary.get("panel_grade", ""),
            "n_barcodes": _coerce_value("n_barcodes", summary.get("n_barcodes", "")),
            "assignment_metric": summary.get("assignment_metric", ""),
            "configured_assignment_k": _coerce_value("configured_assignment_k", summary.get("configured_assignment_k", "")),
            "minimum_hamming_distance": _coerce_value("minimum_hamming_distance", summary.get("minimum_hamming_distance", "")),
            "collision_pairs": _coerce_value("collision_pairs", summary.get("collision_pairs", "")),
            "ambiguous_error_spheres": _coerce_value("ambiguous_error_spheres", summary.get("ambiguous_error_spheres", "")),
            "silent_assignment_risk": _coerce_value("silent_assignment_risk", summary.get("silent_assignment_risk", "")),
            "safe_for_k1_hamming": summary.get("safe_for_k1_hamming", ""),
        }
    }


def parse_top_unmatched_tsv(path: str | Path, limit: int = 10) -> dict[str, dict[str, Any]]:
    """Parse top_unmatched.tsv-style diagnostics, capped for compact reports."""
    path = Path(path)
    data: dict[str, dict[str, Any]] = {}
    for i, row in enumerate(_read_tsv_rows(path)):
        if i >= limit:
            break
        seq = row.get("sequence") or row.get("observed_sequence") or row.get("target") or f"row_{i + 1}"
        data[f"{path.stem}:{seq}"] = row
    return data


def parse_dotmatch_artifacts(paths: list[str | Path]) -> dict[str, dict[str, dict[str, Any]]]:
    """Parse all supported DotMatch report artifacts without importing MultiQC."""
    parsed: dict[str, dict[str, dict[str, Any]]] = {
        "sample_qc": {},
        "summary": {},
        "crispr_qc": {},
        "assay_manifest": {},
        "panel_summary": {},
        "top_unmatched": {},
    }
    for raw_path in paths:
        path = Path(raw_path)
        name = path.name
        if name.endswith("sample_qc.tsv"):
            parsed["sample_qc"].update(parse_sample_qc_tsv(path))
        elif name.endswith("crispr_qc.summary.tsv"):
            parsed["crispr_qc"].update(parse_crispr_qc_summary_tsv(path))
        elif name.endswith("assay_manifest.summary.tsv"):
            parsed["assay_manifest"].update(parse_assay_manifest_summary_tsv(path))
        elif name.endswith("panel_summary.json"):
            parsed["panel_summary"].update(parse_panel_summary_json(path))
        elif name.endswith("top_unmatched.tsv"):
            parsed["top_unmatched"].update(parse_top_unmatched_tsv(path))
        elif name.endswith("summary.json"):
            parsed["summary"].update(parse_summary_json(path))
    return {key: value for key, value in parsed.items() if value}


class DotMatchModule(BaseMultiqcModule):
    """
    DotMatch MultiQC module.

    Add this to a directory and point MultiQC at it, or subclass/inline the
    logic for custom pipelines.
    """

    def __init__(self) -> None:
        if not _HAS_MULTIQC:
            raise ImportError(
                "multiqc package is required to instantiate DotMatchModule. "
                "Install with pip install multiqc and then use the module."
            )
        super().__init__(
            name="DotMatch",
            anchor="dotmatch",
            href="https://github.com/dnncha/dotmatch",
            info="Deterministic short-DNA known-target assignment with explicit ambiguity handling.",
            # Add more metadata as desired
        )

        load_config()

        # Find files
        self.sample_qc_files: list[Any] = list(self.find_log_files(
            "dotmatch/sample_qc", filehandles=False, filecontents=False
        ))
        self.summary_files: list[Any] = list(self.find_log_files(
            "dotmatch/summary_json", filehandles=False, filecontents=False
        ))
        self.crispr_qc_files: list[Any] = list(self.find_log_files(
            "dotmatch/crispr_qc", filehandles=False, filecontents=False
        ))
        self.assay_manifest_files: list[Any] = list(self.find_log_files(
            "dotmatch/assay_manifest", filehandles=False, filecontents=False
        ))
        self.panel_summary_files: list[Any] = list(self.find_log_files(
            "dotmatch/panel_summary", filehandles=False, filecontents=False
        ))
        self.top_unmatched_files: list[Any] = list(self.find_log_files(
            "dotmatch/top_unmatched", filehandles=False, filecontents=False
        ))

        if not any([
            self.sample_qc_files,
            self.summary_files,
            self.crispr_qc_files,
            self.assay_manifest_files,
            self.panel_summary_files,
            self.top_unmatched_files,
        ]):
            raise ModuleNotFoundError("No DotMatch logs found")

        # Parse and add sections
        self._parse_sample_qc()
        self._parse_summary_json()
        self._parse_crispr_qc()
        self._parse_assay_manifest()
        self._parse_panel_summary()
        self._parse_top_unmatched()

        # Cleanup
        for f in (
            self.sample_qc_files
            + self.summary_files
            + self.crispr_qc_files
            + self.assay_manifest_files
            + self.panel_summary_files
            + self.top_unmatched_files
        ):
            self.add_data_source(f)

    def _parse_sample_qc(self) -> None:
        data: dict[str, dict[str, Any]] = {}
        for f in self.sample_qc_files:
            try:
                data.update(parse_sample_qc_tsv(_multiqc_file_path(f)))
            except Exception as exc:
                log.warning(f"Could not parse DotMatch sample_qc {_multiqc_file_path(f)}: {exc}")
                continue

        if data:
            self.add_section(
                name="DotMatch Sample QC",
                anchor="dotmatch-sample-qc",
                description="Per-sample assignment outcomes from DotMatch. Ambiguous reads are excluded from unique counts.",
                plot=_get_table_plot()(data, {
                    "assigned_reads": {"title": "Assigned"},
                    "exact_reads": {"title": "Exact"},
                    "ambiguous_reads": {"title": "Ambiguous"},
                    "no_match_reads": {"title": "No Match"},
                    "assignment_rate": {"title": "Assign Rate"},
                    "ambiguous_rate": {"title": "Ambig Rate"},
                }, {
                    "id": "dotmatch_sample_qc_table",
                    "title": "DotMatch Sample QC",
                }),
            )
            self.general_stats_addcols(data, {
                "assignment_rate": {"title": "Assign %", "description": "Fraction of valid reads uniquely assigned by DotMatch (ambiguous excluded)", "max": 1, "min": 0, "scale": "YlGn"},
                "ambiguous_rate": {"title": "Ambig %", "description": "Fraction of reads that were ambiguous (multiple targets within radius)", "max": 1, "min": 0, "scale": "OrRd"},
            })

    def _parse_summary_json(self) -> None:
        data: dict[str, dict[str, Any]] = {}
        for f in self.summary_files:
            try:
                data.update(parse_summary_json(_multiqc_file_path(f)))
            except Exception as exc:
                log.warning(f"Could not parse DotMatch summary.json {_multiqc_file_path(f)}: {exc}")
                continue

        if data:
            self.add_section(
                name="DotMatch Assignment Summary",
                anchor="dotmatch-assignment-summary",
                description="Run-level assignment outcomes from DotMatch summary.json files.",
                plot=_get_table_plot()(data, {
                    "assigned_unique": {"title": "Unique"},
                    "assigned_exact": {"title": "Exact"},
                    "assigned_corrected": {"title": "Corrected"},
                    "ambiguous": {"title": "Ambiguous"},
                    "unmatched": {"title": "Unmatched"},
                    "invalid": {"title": "Invalid"},
                    "assignment_rate": {"title": "Assign Rate"},
                    "assignment_engine": {"title": "Engine"},
                }, {
                    "id": "dotmatch_assignment_summary_table",
                    "title": "DotMatch Assignment Summary",
                }),
            )

    def _parse_crispr_qc(self) -> None:
        data: dict[str, dict[str, Any]] = {}
        for f in self.crispr_qc_files:
            try:
                data.update(parse_crispr_qc_summary_tsv(_multiqc_file_path(f)))
            except Exception as exc:
                log.warning(f"Could not parse DotMatch crispr_qc.summary {_multiqc_file_path(f)}: {exc}")
                continue

        if data:
            self.add_section(
                name="DotMatch CRISPR QC",
                anchor="dotmatch-crispr-qc",
                description="CRISPR guide-count library QC from DotMatch, including coverage, zero-count guides, Gini, top 1% dominance, and assignment rates.",
                plot=_get_table_plot()(data, {
                    "qc_status": {"title": "QC"},
                    "total_count": {"title": "Total Guides"},
                    "coverage_fraction": {"title": "Coverage"},
                    "zero_count_fraction": {"title": "Zero Frac"},
                    "gini_index": {"title": "Gini"},
                    "top_1pct_fraction": {"title": "Top 1%"},
                    "assignment_rate": {"title": "Assign Rate"},
                }, {
                    "id": "dotmatch_crispr_qc_table",
                    "title": "DotMatch CRISPR QC",
                }),
            )
            self.general_stats_addcols(data, {
                "coverage_fraction": {"title": "Cov %", "description": "Fraction of guide library observed in unique counts", "max": 1, "min": 0, "scale": "YlGn"},
                "zero_count_fraction": {"title": "Zero %", "description": "Fraction of guides with zero unique counts", "max": 1, "min": 0, "scale": "OrRd"},
            })

    def _parse_assay_manifest(self) -> None:
        data: dict[str, dict[str, Any]] = {}
        for f in self.assay_manifest_files:
            try:
                data.update(parse_assay_manifest_summary_tsv(_multiqc_file_path(f)))
            except Exception as exc:
                log.warning(f"Could not parse DotMatch assay_manifest.summary {_multiqc_file_path(f)}: {exc}")
                continue

        if data:
            self.add_section(
                name="DotMatch Assay Manifest",
                anchor="dotmatch-assay-manifest",
                description="Provenance and summary from `dotmatch assay run`, with paths to the full report and manifest.",
                plot=_get_table_plot()(data, {
                    "assay_type": {"title": "Type"},
                    "status": {"title": "Status"},
                    "sample_count": {"title": "Samples"},
                    "autopsy_triggered": {"title": "Autopsy"},
                    "warning_count": {"title": "Warnings"},
                }, {
                    "id": "dotmatch_assay_manifest_table",
                    "title": "DotMatch Assay Manifest",
                }),
            )

    def _parse_top_unmatched(self) -> None:
        data: dict[str, dict[str, Any]] = {}
        for f in self.top_unmatched_files:
            try:
                data.update(parse_top_unmatched_tsv(_multiqc_file_path(f)))
            except Exception as exc:
                log.warning(f"Could not parse DotMatch top_unmatched.tsv {_multiqc_file_path(f)}: {exc}")
                continue

        if data:
            self.add_section(
                name="DotMatch Top Unmatched",
                anchor="dotmatch-top-unmatched",
                description="Most frequent unassigned observed sequences from DotMatch diagnostics.",
                plot=_get_table_plot()(data, None, {
                    "id": "dotmatch_top_unmatched_table",
                    "title": "DotMatch Top Unmatched",
                }),
            )

    def _parse_panel_summary(self) -> None:
        data: dict[str, dict[str, Any]] = {}
        for f in self.panel_summary_files:
            try:
                data.update(parse_panel_summary_json(_multiqc_file_path(f)))
            except Exception as exc:
                log.warning(f"Could not parse DotMatch panel_summary.json {_multiqc_file_path(f)}: {exc}")
                continue

        if data:
            self.add_section(
                name="DotMatch Panel Safety",
                anchor="dotmatch-panel-safety",
                description="Barcode panel check summaries from DotMatch.",
                plot=_get_table_plot()(data, {
                    "status": {"title": "Status"},
                    "panel_grade": {"title": "Grade"},
                    "n_barcodes": {"title": "Barcodes"},
                    "assignment_metric": {"title": "Metric"},
                    "configured_assignment_k": {"title": "k"},
                    "minimum_hamming_distance": {"title": "Min Hamming"},
                    "collision_pairs": {"title": "Collisions"},
                    "ambiguous_error_spheres": {"title": "Ambig Spheres"},
                    "silent_assignment_risk": {"title": "Silent Risk"},
                }, {
                    "id": "dotmatch_panel_safety_table",
                    "title": "DotMatch Panel Safety",
                }),
            )


# For users who want to register without subclassing the whole thing
def load_dotmatch_multiqc() -> None:
    """Helper to make the module discoverable in some MultiQC setups."""
    load_config()
