"""Small, dependency-free sequence and interval primitives."""
from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .models import Edit, Interval

_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


def apply_edits(reference: str, edits: tuple[Edit, ...]) -> str:
    """Apply validated nonoverlapping replacements in original reference coordinates."""
    parts: list[str] = []
    cursor = 0
    for edit in edits:
        parts.extend((reference[cursor:edit.start], edit.sequence))
        cursor = edit.end
    parts.append(reference[cursor:])
    return "".join(parts)


def disrupts_site(edit: Edit, site: Interval) -> bool:
    """Loss of the pristine annotated site, not a thermodynamic prediction of PCR failure.

    Insertion exactly at a site's outer boundary does not disrupt that site.
    Any nonempty replacement overlap is conservatively considered disruptive.
    """
    if edit.start == edit.end:
        return site.start < edit.start < site.end
    return edit.start < site.end and edit.end > site.start


def map_intact_base(position: int, edits: tuple[Edit, ...]) -> int:
    """Map an original surviving base; insertions immediately before it shift it."""
    return position + sum(
        len(edit.sequence) - (edit.end - edit.start)
        for edit in edits if edit.end <= position
    )


def find_all(sequence: str, query: str) -> tuple[int, ...]:
    if not query:
        raise ValueError("empty sequence query")
    hits: list[int] = []
    start = 0
    while (hit := sequence.find(query, start)) != -1:
        hits.append(hit)
        start = hit + 1
    return tuple(hits)
