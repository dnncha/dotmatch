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
- crispr_qc.summary.tsv
- assay_manifest.summary.tsv
and adds nice tables + plots to the MultiQC report.

See examples/workflows/multiqc/ for the simple custom-content alternative
that works today with no extra code.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

log = logging.getLogger(__name__)

# Lazy imports: only fail when actually using the MultiQC integration
# (users who want it do: pip install multiqc "dotmatch[multiqc]" or similar)
# This enables first-class MultiQC support in nf-core / industry pipelines while keeping
# the core package lightweight.
#
# Registration happens via pyproject.toml entry point:
#   [project.entry-points."multiqc.modules.v1"]
#   dotmatch = "dotmatch.multiqc:DotMatchModule"
# After `pip install "dotmatch[multiqc]"`, `multiqc ...` should auto-discover it
# (or use --module dotmatch explicitly).
_HAS_MULTIQC = False
try:
    from multiqc import BaseMultiqcModule  # type: ignore
    from multiqc.plots import table as _mqc_table  # type: ignore
    _HAS_MULTIQC = True
except Exception:  # noqa: BLE001
    BaseMultiqcModule = object  # type: ignore
    _mqc_table = None  # type: ignore

def _get_table_plot():
    if not _HAS_MULTIQC or _mqc_table is None:
        raise ImportError("multiqc is required to use the DotMatch MultiQC module")
    return _mqc_table.plot


# Pure, dependency-free parsers. These can be used from any Python code (notebooks,
# custom pipelines, other report generators) for excellent UX and to guarantee
# that the exact DotMatch semantics (unique-only counts, documented rates, etc.)
# are preserved when consuming the TSV artifacts outside MultiQC.
def parse_sample_qc_tsv(path: str | Path) -> dict[str, dict[str, Any]]:
    """Parse DotMatch sample_qc.tsv into dict of per-sample metrics.

    Columns and semantics exactly as documented in docs/schemas.md.
    assignment_rate etc. denominator is valid_extracted_reads; ambiguous reads
    are never counted toward any target (core scientific accuracy property).
    """
    path = Path(path)
    data: dict[str, dict[str, Any]] = {}
    with path.open(encoding="utf-8") as fh:
        header = fh.readline().strip().split("\t")
        for line in fh:
            if not line.strip():
                continue
            vals = line.strip().split("\t")
            row = dict(zip(header, vals))
            sname = row.get("sample_id", row.get("sample", path.stem))
            data[sname] = {
                "total_reads": row.get("total_reads", ""),
                "assigned_reads": row.get("assigned_reads", ""),
                "exact_reads": row.get("exact_reads", ""),
                "ambiguous_reads": row.get("ambiguous_reads", ""),
                "no_match_reads": row.get("no_match_reads", ""),
                "invalid_reads": row.get("invalid_reads", ""),
                "assignment_rate": row.get("assignment_rate", ""),
                "ambiguous_rate": row.get("ambiguous_rate", ""),
                "gini_index": row.get("gini_index", ""),
                "candidates_verified": row.get("candidates_verified", ""),
            }
    return data

def parse_crispr_qc_summary_tsv(path: str | Path) -> dict[str, dict[str, Any]]:
    """Parse crispr_qc.summary.tsv (guide library representation + QC)."""
    path = Path(path)
    data: dict[str, dict[str, Any]] = {}
    with path.open(encoding="utf-8") as fh:
        header = fh.readline().strip().split("\t")
        for line in fh:
            if not line.strip():
                continue
            vals = line.strip().split("\t")
            row = dict(zip(header, vals))
            sname = row.get("sample_id", path.stem)
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
    with path.open(encoding="utf-8") as fh:
        header = fh.readline().strip().split("\t")
        for line in fh:
            if not line.strip():
                continue
            vals = line.strip().split("\t")
            row = dict(zip(header, vals))
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

        # Find files
        self.sample_qc_files: list[Path] = self.find_log_files(
            "dotmatch/sample_qc", filehandles=False, filecontents=False
        )
        self.crispr_qc_files: list[Path] = self.find_log_files(
            "dotmatch/crispr_qc", filehandles=False, filecontents=False
        )
        self.assay_manifest_files: list[Path] = self.find_log_files(
            "dotmatch/assay_manifest", filehandles=False, filecontents=False
        )

        if not any([self.sample_qc_files, self.crispr_qc_files, self.assay_manifest_files]):
            raise ModuleNotFoundError("No DotMatch logs found")

        # Parse and add sections
        self._parse_sample_qc()
        self._parse_crispr_qc()
        self._parse_assay_manifest()

        # Cleanup
        for f in self.sample_qc_files + self.crispr_qc_files + self.assay_manifest_files:
            self.add_data_source(f)

    def _parse_sample_qc(self) -> None:
        data: dict[str, dict[str, Any]] = {}
        for f in self.sample_qc_files:
            try:
                data.update(parse_sample_qc_tsv(f["fn"]))
            except Exception as exc:
                log.warning(f"Could not parse DotMatch sample_qc {f['fn']}: {exc}")
                continue

        if data:
            self.add_section(
                name="DotMatch Sample QC",
                anchor="dotmatch-sample-qc",
                description="Per-sample assignment outcomes from DotMatch. Only uniquely assigned reads contribute to counts (ambiguous reads are excluded by design for scientific accuracy).",
                plot=_get_table_plot()(data, {
                    "assigned_reads": {"title": "Assigned"},
                    "exact_reads": {"title": "Exact"},
                    "ambiguous_reads": {"title": "Ambiguous"},
                    "no_match_reads": {"title": "No Match"},
                    "assignment_rate": {"title": "Assign Rate"},
                    "ambiguous_rate": {"title": "Ambig Rate"},
                }),
            )
            self.general_stats_addcols(data, {
                "assignment_rate": {"title": "Assign %", "description": "Fraction of valid reads uniquely assigned by DotMatch (ambiguous excluded)", "max": 1, "min": 0, "scale": "YlGn"},
                "ambiguous_rate": {"title": "Ambig %", "description": "Fraction of reads that were ambiguous (multiple targets within radius)", "max": 1, "min": 0, "scale": "OrRd"},
            })

    def _parse_crispr_qc(self) -> None:
        data: dict[str, dict[str, Any]] = {}
        for f in self.crispr_qc_files:
            try:
                data.update(parse_crispr_qc_summary_tsv(f["fn"]))
            except Exception as exc:
                log.warning(f"Could not parse DotMatch crispr_qc.summary {f['fn']}: {exc}")
                continue

        if data:
            self.add_section(
                name="DotMatch CRISPR QC",
                anchor="dotmatch-crispr-qc",
                description="CRISPR guide-count library QC and representation from DotMatch (zero-count guides, Gini, top 1% dominance, assignment rates). High zero-count or high Gini indicates poor library coverage or skew — important for screen interpretability.",
                plot=_get_table_plot()(data, {
                    "qc_status": {"title": "QC"},
                    "total_count": {"title": "Total Guides"},
                    "coverage_fraction": {"title": "Coverage"},
                    "zero_count_fraction": {"title": "Zero Frac"},
                    "gini_index": {"title": "Gini"},
                    "top_1pct_fraction": {"title": "Top 1%"},
                    "assignment_rate": {"title": "Assign Rate"},
                }),
            )
            self.general_stats_addcols(data, {
                "coverage_fraction": {"title": "Cov %", "description": "Fraction of guide library observed in unique counts", "max": 1, "min": 0, "scale": "YlGn"},
                "zero_count_fraction": {"title": "Zero %", "description": "Fraction of guides with zero unique counts (bad for screens)", "max": 1, "min": 0, "scale": "OrRd"},
            })

    def _parse_assay_manifest(self) -> None:
        data: dict[str, dict[str, Any]] = {}
        for f in self.assay_manifest_files:
            try:
                data.update(parse_assay_manifest_summary_tsv(f["fn"]))
            except Exception as exc:
                log.warning(f"Could not parse DotMatch assay_manifest.summary {f['fn']}: {exc}")
                continue

        if data:
            self.add_section(
                name="DotMatch Assay Manifest",
                anchor="dotmatch-assay-manifest",
                description="Provenance and summary from `dotmatch assay run`. Links to full HTML report and JSON for full audit trail and reproducibility (critical for scientific accuracy and industry use).",
                plot=_get_table_plot()(data, {
                    "assay_type": {"title": "Type"},
                    "status": {"title": "Status"},
                    "sample_count": {"title": "Samples"},
                    "autopsy_triggered": {"title": "Autopsy"},
                    "warning_count": {"title": "Warnings"},
                }),
            )


# For users who want to register without subclassing the whole thing
def load_dotmatch_multiqc() -> None:
    """Helper to make the module discoverable in some MultiQC setups."""
    # In practice, MultiQC discovers via package entry points or --module path
    pass
