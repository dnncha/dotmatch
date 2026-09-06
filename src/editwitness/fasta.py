"""Local single-record FASTA import and unambiguous exact primer-site location."""
from __future__ import annotations

from pathlib import Path

from ._version import EXACT_MODEL_VERSION
from .io import MAX_INPUT_BYTES, InputError
from .models import Allele, Assay, Edit, Hypothesis, Interval, Manifest, Reference
from .sequence import find_all, reverse_complement


def read_fasta(path: str | Path) -> tuple[str, str]:
    with Path(path).open("rb") as handle:
        data = handle.read(MAX_INPUT_BYTES + 1)
    if len(data) > MAX_INPUT_BYTES:
        raise InputError("FASTA exceeds input byte budget")
    try:
        lines = [line.strip() for line in data.decode("utf-8").splitlines() if line.strip()]
    except UnicodeDecodeError as error:
        raise InputError("FASTA must be UTF-8 text") from error
    if not lines or not lines[0].startswith(">") or not lines[0][1:].strip():
        raise InputError("FASTA must have a nonempty header")
    if any(line.startswith(">") for line in lines[1:]):
        raise InputError("exactly one local-reference FASTA record is required")
    sequence = "".join(lines[1:]).upper()
    if not sequence or set(sequence) - set("ACGT"):
        raise InputError("reference must contain only A/C/G/T; ambiguity bases are not modeled")
    if len(sequence) > 20_000:
        raise InputError("use a local reference window of at most 20,000 bases")
    return lines[0][1:], sequence


def init_from_fasta(path: str | Path, left_oligo: str, right_oligo: str,
                    edit_position: int, alternate: str, *, read_bases: int | None = None) -> Manifest:
    name, sequence = read_fasta(path)
    left_oligo, right_oligo, alternate = left_oligo.upper(), right_oligo.upper(), alternate.upper()
    for oligo in (left_oligo, right_oligo):
        if not oligo or set(oligo) - set("ACGT"):
            raise InputError("primers must contain only A/C/G/T")
    left_hits = find_all(sequence, left_oligo)
    right_hits = find_all(sequence, reverse_complement(right_oligo))
    if len(left_hits) != 1 or len(right_hits) != 1:
        raise InputError("each primer must have exactly one inward-orientation exact local match; "
                         "use explicit annotated coordinates for a deliberate alternative")
    if len(alternate) != 1 or alternate not in "ACGT" or not 0 <= edit_position < len(sequence):
        raise InputError("init requires one alternate A/C/G/T base and a valid local 0-based edit position")
    left = Interval(start=left_hits[0], end=left_hits[0] + len(left_oligo))
    right = Interval(start=right_hits[0], end=right_hits[0] + len(right_oligo))
    if not left.end <= edit_position < right.start:
        raise InputError("init requires the intended substitution between the primer binding sites")
    return Manifest(
        observation_model=EXACT_MODEL_VERSION,
        reference=Reference(name=name, sequence=sequence),
        alleles=(Allele(id="reference"), Allele(id="intended", edits=(
            Edit(start=edit_position, end=edit_position + 1, sequence=alternate),
        )), Allele(id="window_deleted", description="Deliberate blind-spot hypothesis, not an observed event",
                  edits=(Edit(start=0, end=len(sequence)),))),
        hypotheses=(Hypothesis(id="intended_biallelic", alleles=("intended", "intended")),
                    Hypothesis(id="intended_reference", alleles=("intended", "reference")),
                    Hypothesis(id="intended_window_deleted", alleles=("intended", "window_deleted"))),
        expected_hypothesis="intended_biallelic",
        assays=(Assay(id="amplicon", left_primer=left, right_primer=right,
                      left_oligo=left_oligo, right_oligo=right_oligo,
                      readout="full_insert" if read_bases is None else "paired_end", read_bases=read_bases),),
    )
