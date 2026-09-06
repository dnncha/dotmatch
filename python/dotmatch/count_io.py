"""Lossless raw-count tables. Missing values are errors, never zero counts."""
from __future__ import annotations

import csv
import gzip
import re
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Sequence

from ._validation import identifier

ID_COLUMNS = ("sgrna", "target_id", "guide_id", "feature_id", "guide", "id", "name")
METADATA_COLUMNS = frozenset({"target_seq", "sequence", "seq", "gene", "gene_id", "gene_symbol", "ambiguous_nearby"})


def parse_count_value(text: str, guide_id: str, sample: str) -> int:
    """Parse an exact nonnegative integer, including integral decimal notation.

    Decimal parsing avoids float64 rounding above 2**53. Limit exponent/input
    growth before integer conversion; even 1,000-digit counts are far beyond
    sequencing workloads, while unbounded exponents can exhaust memory.
    """
    location = f"{guide_id}/{sample}"
    if not isinstance(text, str) or not text.strip():
        raise ValueError(f"missing count for {location}; use an explicit 0")
    text = text.strip()
    if len(text) > 4096:
        raise ValueError(f"count too large for {location}")
    if re.fullmatch(r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?", text) is None:
        raise ValueError(f"non-numeric or non-finite count for {location}")
    try:
        value = Decimal(text)
    except InvalidOperation as exc:
        raise ValueError(f"non-numeric count for {location}") from exc
    if not value.is_finite():
        raise ValueError(f"non-finite count for {location}")
    if value < 0:
        raise ValueError(f"negative count for {location}")
    if value != value.to_integral_value():
        raise ValueError(f"non-integer count for {location}")
    if value.adjusted() >= 1000:
        raise ValueError(f"count too large for {location}; at most 1,000 digits are supported")
    return int(value)


@dataclass(frozen=True)
class CountTable:
    id_column: str
    target_ids: tuple[str, ...]
    sample_columns: tuple[str, ...]
    sample_names: tuple[str, ...]
    counts: tuple[tuple[int, ...], ...]  # targets x samples, raw integer counts
    metadata: dict[str, tuple[str, ...]]


def read_count_table(path: str | Path, *, sample_cols: Sequence[str] | None = None,
                     var_cols: Sequence[str] = (), mageck_only: bool = False) -> CountTable:
    """Read MAGeCK or DotMatch TSV without guessing numeric columns from values.

    DotMatch detailed output contributes only *_count_total columns by default;
    exact/corrected components and QC fields must not become additional samples.
    sample_cols names source columns and preserves the requested order.
    """
    source = Path(path)
    opener = gzip.open if source.name.lower().endswith(".gz") else open
    with opener(source, "rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t", strict=True)
        header = next(reader, None)
        if header is None or len(header) < (3 if mageck_only else 2):
            raise ValueError("count matrix must have guide, gene, and at least one sample column" if mageck_only else "count table needs an ID and at least one count column")
        if any(not name or name != name.strip() or any(ord(c) < 32 or ord(c) == 127 for c in name) for name in header):
            raise ValueError("count matrix columns must be nonempty without surrounding whitespace or controls")
        if len(set(header)) != len(header) or len({name.casefold() for name in header}) != len(header):
            raise ValueError("count matrix columns must be unique")
        lower = {name.lower(): name for name in header}
        id_col = header[0] if mageck_only else next((lower[name] for name in ID_COLUMNS if name in lower), header[0])
        metadata_cols = [name for name in header if name != id_col and (name.lower() in METADATA_COLUMNS or name in var_cols)]
        detailed = not mageck_only and "target_seq" in lower and any(name.endswith("_count_total") for name in header)
        if mageck_only:
            metadata_cols = [header[1]]
            candidates = header[2:]
        elif detailed:
            candidates = [name for name in header if name.endswith("_count_total")]
        else:
            candidates = [name for name in header if name != id_col and name not in metadata_cols]
        selected = list(candidates) if sample_cols is None else list(sample_cols)
        if not selected or len(set(selected)) != len(selected) or any(name not in candidates for name in selected):
            raise ValueError("sample_cols must contain distinct existing count-column names")
        sample_names = [name[:-len("_count_total")] if detailed else name for name in selected]
        if any(not name for name in sample_names) or len(set(sample_names)) != len(sample_names):
            raise ValueError("sample names derived from count columns must be nonempty and unique")
        ids, counts, metadata, seen = [], [], {name: [] for name in metadata_cols}, set()
        for fields in reader:
            if len(fields) != len(header):
                raise ValueError(f"count matrix row {reader.line_num} has {len(fields)} fields; expected {len(header)}")
            row = dict(zip(header, fields))
            target_id = identifier(row[id_col], "guide id")
            if target_id in seen:
                raise ValueError(f"count matrix contains duplicate guide id: {target_id}")
            seen.add(target_id)
            # Validate all candidate columns, even when the caller selects a subset.
            values = {name: parse_count_value(row[name], target_id, name) for name in candidates}
            ids.append(target_id)
            counts.append(tuple(values[name] for name in selected))
            for name in metadata:
                value = row[name]
                if any(ord(ch) < 32 or ord(ch) == 127 for ch in value):
                    raise ValueError(f"control character in {name} for {target_id}")
                metadata[name].append(value)
        if not ids:
            raise ValueError("count matrix contains no guides")
    return CountTable(id_col, tuple(ids), tuple(selected), tuple(sample_names), tuple(counts),
                      {key: tuple(values) for key, values in metadata.items()})


def read_crispr_count_matrix(path: str | Path) -> dict[str, object]:
    """Preserve the established CLI/QC dictionary contract."""
    table = read_count_table(path, mageck_only=True)
    gene_col = next(iter(table.metadata))
    guides = [{"id": target_id, "gene": table.metadata[gene_col][i],
               "counts": dict(zip(table.sample_names, table.counts[i]))}
              for i, target_id in enumerate(table.target_ids)]
    return {"guide_col": table.id_column, "gene_col": gene_col,
            "samples": list(table.sample_names), "guides": guides,
            "sample_counts": {name: [row[i] for row in table.counts] for i, name in enumerate(table.sample_names)}}
