"""Versioned, strict, immutable data contracts. Coordinates are local and half-open."""
from __future__ import annotations

from typing import Annotated, Literal, Self, TypeVar

from pydantic import BaseModel, BeforeValidator, ConfigDict, Field, model_validator

T = TypeVar("T")


def _tuple(value: object) -> object:
    return tuple(value) if isinstance(value, list) else value


Items = Annotated[tuple[T, ...], BeforeValidator(_tuple)]
Identifier = Annotated[str, Field(min_length=1, max_length=80, pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")]
DNA = Annotated[str, Field(pattern=r"^[ACGT]*$", max_length=20_000)]


class Contract(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True, strict=True, allow_inf_nan=False, revalidate_instances="always")


class Interval(Contract):
    """Local reference interval [start, end), not a VCF or genomic position."""
    start: int = Field(ge=0)
    end: int = Field(gt=0)

    @model_validator(mode="after")
    def ordered(self) -> Self:
        if self.start >= self.end:
            raise ValueError("interval requires start < end")
        return self


class Reference(Contract):
    name: str = Field(min_length=1, max_length=200)
    sequence: DNA = Field(min_length=1)
    assembly: str | None = Field(default=None, max_length=200)
    contig: str | None = Field(default=None, max_length=200)
    genomic_start: int | None = Field(default=None, ge=0)
    synthetic: bool = False


class Edit(Contract):
    """Replace reference[start:end] with sequence; start == end is an insertion."""
    start: int = Field(ge=0)
    end: int = Field(ge=0)
    sequence: DNA = ""

    @model_validator(mode="after")
    def ordered(self) -> Self:
        if self.end < self.start:
            raise ValueError("edit requires start <= end")
        if self.start == self.end and not self.sequence:
            raise ValueError("empty insertion is not an edit")
        return self


class Allele(Contract):
    id: Identifier
    description: str = Field(default="", max_length=1000)
    edits: Items[Edit] = Field(default=(), max_length=64)


class Assay(Contract):
    id: Identifier
    description: str = Field(default="", max_length=1000)
    left_primer: Interval
    right_primer: Interval
    left_oligo: DNA | None = Field(default=None, min_length=1)
    right_oligo: DNA | None = Field(default=None, min_length=1)
    readout: Literal["full_insert", "paired_end"] = "full_insert"
    read_bases: int | None = Field(default=None, ge=1, le=20_000)
    min_product_bp: int = Field(default=1, ge=1)
    max_product_bp: int | None = Field(default=None, ge=1)
    cost_units: int = Field(default=1, ge=1, le=1_000_000)

    @model_validator(mode="after")
    def consistent(self) -> Self:
        if self.left_primer.end >= self.right_primer.start:
            raise ValueError("primers must face inward with nonempty intervening reference sequence")
        if (self.readout == "paired_end") != (self.read_bases is not None):
            raise ValueError("read_bases is required only for paired_end; it counts post-trim insert bases")
        if self.max_product_bp is not None and self.max_product_bp < self.min_product_bp:
            raise ValueError("max_product_bp must be >= min_product_bp")
        return self


class Hypothesis(Contract):
    id: Identifier
    description: str = Field(default="", max_length=1000)
    alleles: Items[Identifier] = Field(min_length=2, max_length=2)


class DeletionScan(Contract):
    """Inclusive endpoint ranges, sampled with step, describing deletions [start,end)."""
    start_min: int = Field(ge=0)
    start_max: int = Field(ge=0)
    end_min: int = Field(ge=1)
    end_max: int = Field(ge=1)
    step: int = Field(default=1, ge=1)
    min_length: int = Field(default=1, ge=1)
    max_length: int | None = Field(default=None, ge=1)

    @model_validator(mode="after")
    def consistent(self) -> Self:
        if self.start_min > self.start_max or self.end_min > self.end_max:
            raise ValueError("scan endpoint minima must not exceed maxima")
        if self.max_length is not None and self.max_length < self.min_length:
            raise ValueError("scan max_length must be >= min_length")
        cells = ((self.start_max - self.start_min) // self.step + 1) * (
            (self.end_max - self.end_min) // self.step + 1
        )
        if cells > 500_000:
            raise ValueError("scan grid exceeds 500,000 endpoint pairs; increase step or narrow ranges")
        return self


class DeletionGeneration(Contract):
    method: Literal["single-reference-deletions-v1"] = "single-reference-deletions-v1"
    input_manifest_sha256: str
    grid: DeletionScan
    fixed_allele: Identifier
    valid_deletions: int = Field(ge=1)
    deduplicated_states: int = Field(ge=0)
    added_alleles: int = Field(ge=0)
    added_hypotheses: int = Field(ge=0)
    max_new_hypotheses: int = Field(ge=1, le=999)
    caveat: str = ("One expected allele fixed; single deletions on the reference haplotype only. "
                   "Equivalent local sequence states deduplicated; first geometry retained. "
                   "No probabilities, sampling or exhaustive biological coverage.")


class Manifest(Contract):
    schema_version: Literal["1.0", "1.1"] = "1.1"
    observation_model: Literal[
        "original-sites-presence-v1", "exact-local-sequence-presence-v2"
    ] = "original-sites-presence-v1"
    coordinate_system: Literal["0-based-half-open"] = "0-based-half-open"
    reference: Reference
    alleles: Items[Allele] = Field(min_length=1, max_length=128)
    hypotheses: Items[Hypothesis] = Field(min_length=1, max_length=1000)
    expected_hypothesis: Identifier
    assays: Items[Assay] = Field(min_length=1, max_length=16)
    candidates: Items[Assay] = Field(default=(), max_length=24)
    deletion_scan: DeletionScan | None = None
    generation: DeletionGeneration | None = None

    @model_validator(mode="after")
    def cross_validate(self) -> Self:
        from .sequence import reverse_complement

        if self.schema_version == "1.0" and self.observation_model != "original-sites-presence-v1":
            raise ValueError("exact-local model requires manifest schema_version 1.1")
        n = len(self.reference.sequence)
        for name, values in (
            ("allele", self.alleles), ("hypothesis", self.hypotheses),
            ("assay", self.assays + self.candidates),
        ):
            ids = [item.id for item in values]
            if len(set(ids)) != len(ids):
                raise ValueError(f"duplicate {name} id")
        known = {a.id for a in self.alleles}
        for allele in self.alleles:
            previous: Edit | None = None
            for edit in allele.edits:
                if edit.end > n:
                    raise ValueError(f"{allele.id}: edit outside reference")
                if self.reference.sequence[edit.start:edit.end] == edit.sequence:
                    raise ValueError(f"{allele.id}: no-op edit; omit it")
                if previous is not None:
                    if edit.start < previous.end or edit.start == previous.start:
                        raise ValueError(f"{allele.id}: edits must be sorted, nonoverlapping and unambiguous")
                previous = edit
        for hypothesis in self.hypotheses:
            if not set(hypothesis.alleles) <= known:
                raise ValueError(f"{hypothesis.id}: unknown allele id")
        if self.expected_hypothesis not in {h.id for h in self.hypotheses}:
            raise ValueError("expected_hypothesis does not name a declared hypothesis")
        for assay in self.assays + self.candidates:
            if assay.right_primer.end > n:
                raise ValueError(f"{assay.id}: primer outside reference")
            left = self.reference.sequence[assay.left_primer.start:assay.left_primer.end]
            right = reverse_complement(
                self.reference.sequence[assay.right_primer.start:assay.right_primer.end]
            )
            if self.observation_model == "exact-local-sequence-presence-v2" and left == right:
                raise ValueError(f"{assay.id}: identical primer oligos have ambiguous read orientation; unsupported")
            if assay.left_oligo is not None and assay.left_oligo != left:
                raise ValueError(f"{assay.id}: left_oligo does not match its annotated reference site")
            if assay.right_oligo is not None and assay.right_oligo != right:
                raise ValueError(f"{assay.id}: right_oligo must be the 5'-to-3' reverse complement")
        reconstructed_bases = sum(n + sum(len(e.sequence) - (e.end - e.start) for e in a.edits)
                                  for a in self.alleles)
        if max(n * len(self.alleles), reconstructed_bases) * (len(self.assays) + len(self.candidates)) > 20_000_000:
            raise ValueError("analysis exceeds conservative 20-million-base work budget; split the manifest")
        if self.deletion_scan is not None:
            if self.deletion_scan.start_max >= n or self.deletion_scan.end_max > n:
                raise ValueError("scan endpoints outside reference")
        return self


class Notice(Contract):
    code: str
    message: str
    related_ids: Items[str] = ()


class ProductObservation(Contract):
    """One exact local heteroprimer product. Coordinates refer to the final allele."""
    plus_left_site: Interval
    plus_right_site: Interval
    orientation: Literal["forward", "reverse"]
    product_length: int = Field(ge=1)
    reads: Items[str]
    signal_id: str


class AlleleEvidence(Contract):
    """Inspectable edit definition plus final-sequence identity, without inferred dosage."""
    allele_id: str
    description: str
    edits: Items[Edit]
    sequence_length: int = Field(ge=0)
    sequence_sha256: str


class AlleleObservation(Contract):
    allele_id: str
    assay_id: str
    status: Literal["potentially_observable", "original_binding_site_disrupted", "outside_product_bounds", "no_exact_local_product"]
    reason: str
    product_length: int | None = None
    reads: Items[str] = ()
    signal_id: str | None = None
    # Singular fields above are populated only for an unambiguous signal; use these collections.
    signal_ids: Items[str] = ()
    products: Items[ProductObservation] = ()


class HypothesisObservation(Contract):
    hypothesis_id: str
    assay_id: str
    signal_ids: Items[str]


class HypothesisAssessment(Contract):
    hypothesis_id: str
    alleles: Items[str]
    equivalent_to_expected: bool
    distinguishing_existing_assays: Items[str]
    same_local_genomic_state_as_expected: bool = False
    description: str = ""


class Witness(Contract):
    hypothesis_id: str
    expected_alleles: Items[str]
    alternative_alleles: Items[str]
    explanation: str
    resolving_candidate_assays: Items[str]


class PanelPlan(Contract):
    goal: Literal["distinguish_expected_from_currently_equivalent_alternatives"] = (
        "distinguish_expected_from_currently_equivalent_alternatives"
    )
    algorithm: Literal["exhaustive_minimum_cost", "greedy_weighted_cover", "not_needed", "no_separating_candidates"]
    optimality: Literal["proven_within_declared_candidates", "not_proven", "not_applicable"]
    selected_assays: Items[str]
    cost_units: int
    resolved_hypotheses: Items[str]
    unresolved_hypotheses: Items[str]
    note: str
    dominated_candidates: Items[str] = ()


class ReferenceSummary(Contract):
    name: str
    length: int
    sequence_sha256: str
    synthetic: bool
    assembly: str | None
    contig: str | None
    genomic_start: int | None


class Analysis(Contract):
    schema_version: Literal["1.0", "1.1"] = "1.1"
    kind: Literal["editwitness.analysis"] = "editwitness.analysis"
    package_version: str
    model_version: str
    manifest_sha256: str
    result_sha256: str = ""
    validation_status: Literal["software-tested; not empirically validated"] = (
        "software-tested; not empirically validated"
    )
    conclusion: Literal["ambiguity_demonstrated", "distinguishable_only_within_declared_model"]
    expected_hypothesis: str
    expected_alleles: Items[str]
    reference: ReferenceSummary
    assays: Items[Assay]
    candidates: Items[Assay]
    assumptions: Items[str]
    notices: Items[Notice]
    generation: DeletionGeneration | None = None
    allele_evidence: Items[AlleleEvidence] = ()
    allele_observations: Items[AlleleObservation]
    hypothesis_observations: Items[HypothesisObservation]
    hypotheses: Items[HypothesisAssessment]
    witnesses: Items[Witness]
    plan: PanelPlan


class ScanAssayCounts(Contract):
    assay_id: str
    potentially_amplifiable: int
    binding_site_disrupted: int
    outside_product_bounds: int


class ScanExample(Contract):
    start: int
    end: int
    length: int
    assay_statuses: dict[str, str]


class ScanResult(Contract):
    schema_version: Literal["1.0", "1.1"] = "1.1"
    kind: Literal["editwitness.deletion_scan"] = "editwitness.deletion_scan"
    package_version: str
    model_version: str
    manifest_sha256: str
    result_sha256: str = ""
    grid: DeletionScan
    enumerated_deletions: int
    all_existing_assays_structurally_blind: int
    assays: Items[ScanAssayCounts]
    blind_examples: Items[ScanExample]
    caveat: str
