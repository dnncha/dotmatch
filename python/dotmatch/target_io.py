"""Target-table input shared by the Python API and Python CLI workflows."""

from __future__ import annotations

import csv
import gzip
from dataclasses import dataclass
from pathlib import Path

ID_COLUMNS = ("target_id", "guide_id", "barcode_id", "id", "name", "sgrna", "guide", "sgrnaid", "sgrna_id")
SEQ_COLUMNS = ("target_seq", "guide_seq", "barcode_seq", "sequence", "seq", "grna.sequence", "bases", "sgrna.sequence", "sgrna_sequence", "sgrna_seq", "guide_sequence", "guidesequence")
GENE_COLUMNS = ("gene", "gene_id", "gene_symbol", "gene.symbol", "target_gene")


@dataclass(frozen=True)
class TargetRecord:
    target_id: str
    sequence: str
    gene: str = ""


def read_target_table(path: str | Path) -> list[TargetRecord]:
    """Read plain/gzipped CSV or TSV, retaining row order and duplicate sequences.

    Named columns may be reordered. IDs must be nonempty and unique: silently
    merging counts under duplicate IDs would conceal an input error. Headerless
    one-column tables retain the historical target_0, target_1, ... IDs. N/IUPAC
    symbols remain literal symbols; no expansion or sequence repair is applied.
    """
    source = Path(path)
    name = source.name.lower()
    compressed = name.endswith(".gz")
    if compressed:
        name = name[:-3]
    delimiter = "," if name.endswith(".csv") else "\t"
    opener = gzip.open if compressed else open
    targets: list[TargetRecord] = []
    seen: set[str] = set()
    first = True
    named = False
    id_col, seq_col, gene_col = 0, 1, 2
    sequence_only = False
    width = None
    with opener(source, "rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter=delimiter, strict=True)
        for row in reader:
            cols = [cell.strip() for cell in row]
            if not cols or not any(cols) or cols[0].startswith("#"):
                continue
            if width is None:
                width = len(cols)
            if len(cols) != width:
                raise ValueError(f"row has a different number of fields in {source} at line {reader.line_num}")
            if first:
                header = [cell.lower() for cell in cols]
                sequence_only = len(header) == 1 and header[0] in SEQ_COLUMNS
                named = sequence_only or (
                    any(c in header for c in ID_COLUMNS)
                    and any(c in header for c in SEQ_COLUMNS)
                )
                first = False
                if named:
                    if any(sum(name in group for name in header) > 1 for group in (ID_COLUMNS, SEQ_COLUMNS, GENE_COLUMNS)):
                        raise ValueError(f"multiple possible ID, sequence or gene columns in {source}; provide an unambiguous library")
                    if len(set(header)) != len(header) or any(not name for name in header):
                        raise ValueError(f"duplicate header columns in {source}")
                    if sequence_only:
                        seq_col = 0
                    else:
                        id_col = next(
                            header.index(c) for c in ID_COLUMNS if c in header
                        )
                        seq_col = next(
                            header.index(c) for c in SEQ_COLUMNS if c in header
                        )
                    gene_col = next(
                        (header.index(c) for c in GENE_COLUMNS if c in header), -1
                    )
                    continue
            if sequence_only or (not named and len(cols) == 1):
                if len(cols) != 1:
                    raise ValueError(
                        f"expected one sequence column in {source} at line {reader.line_num}"
                    )
                target_id, sequence, gene = f"target_{len(targets)}", cols[0], ""
            else:
                required = max(id_col, seq_col, gene_col if named else -1)
                if required >= len(cols):
                    raise ValueError(
                        f"missing target columns in {source} at line {reader.line_num}"
                    )
                target_id, sequence = cols[id_col], cols[seq_col]
                gene = cols[gene_col] if 0 <= gene_col < len(cols) else ""
            if target_id.startswith('"') or gene.startswith('"'):
                raise ValueError("target ID and gene may not start with a literal double quote (TSV identity boundary)")
            if not target_id:
                raise ValueError(
                    f"empty target ID in {source} at line {reader.line_num}"
                )
            if target_id in seen:
                raise ValueError(
                    f"duplicate target ID {target_id!r} in {source}; use distinct IDs, even for duplicate sequences"
                )
            if not sequence:
                raise ValueError(
                    f"empty target sequence in {source} at line {reader.line_num}"
                )
            if not sequence.isascii() or any(
                ch.isspace() or ord(ch) < 32 or ord(ch) == 127 for ch in sequence
            ):
                raise ValueError(
                    f"target sequence must contain ASCII symbols without whitespace in {source} at line {reader.line_num}"
                )
            if any(ord(ch) < 32 or ord(ch) == 127 for ch in target_id + gene):
                raise ValueError(
                    f"target ID and gene must not contain control characters in {source} at line {reader.line_num}"
                )
            seen.add(target_id)
            targets.append(TargetRecord(target_id, sequence.upper(), gene))
    if not targets:
        raise ValueError(f"no targets found in {source}")
    return targets
