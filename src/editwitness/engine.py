"""Pure, deterministic analysis; the report layer never changes scientific conclusions."""
from __future__ import annotations

import hashlib

from ._version import MODEL_VERSION, __version__
from .io import digest, seal
from .models import (
    Analysis, HypothesisAssessment, HypothesisObservation, Manifest, Notice,
    ReferenceSummary, Witness,
)
from .observations import observe_allele
from .planner import plan_panel
from .sequence import find_all

ASSUMPTIONS = (
    "Only the explicitly declared local alleles and diploid clonal hypotheses are compared. "
    "No exhaustive biological-outcome coverage is claimed.",
    "Amplification uses only the annotated original primer sites, which must remain pristine. "
    "Primer mismatches, rescued or newly formed binding sites and nonspecific amplification are not modeled.",
    "Every eligible product is detected without error. Real sampling, PCR bias, assay failure and "
    "limits of detection are not modeled.",
    "An assay observes the SET of full-insert sequences or ordered post-primer-trim read pairs. "
    "Allele multiplicity, read counts and genomic copy number are deliberately not observed.",
    "Configured product-size bounds are hard user-declared inclusion bounds, not a model of PCR efficiency.",
    "Paired-end readout does not observe the unsequenced gap or insert length. Full-insert mode "
    "requires genuinely observing the entire primer-trimmed insert, not merely using a long-read platform.",
    "Model counterexamples are not evidence that a defect actually occurred. No clinical or clone-release "
    "decision is authorized by this software.",
)


def analyze(manifest: Manifest) -> Analysis:
    """Compare declared hypotheses and return explicit observational counterexamples.

    The input must be a validated Manifest (use load_manifest for files). This function
    has no side effects, performs no network access, and never modifies caller outputs.
    """
    assays = manifest.assays + manifest.candidates
    expected = next(h for h in manifest.hypotheses if h.id == manifest.expected_hypothesis)
    observations = tuple(
        observe_allele(manifest.reference.sequence, allele, assay)
        for allele in manifest.alleles for assay in assays
    )
    by_pair = {(o.allele_id, o.assay_id): o for o in observations}
    h_observations: list[HypothesisObservation] = []
    signatures: dict[tuple[str, str], tuple[str, ...]] = {}
    for hypothesis in manifest.hypotheses:
        for assay in assays:
            ids = tuple(sorted({
                signal for allele_id in hypothesis.alleles
                if (signal := by_pair[allele_id, assay.id].signal_id) is not None
            }))
            signatures[hypothesis.id, assay.id] = ids
            h_observations.append(HypothesisObservation(
                hypothesis_id=hypothesis.id, assay_id=assay.id, signal_ids=ids
            ))
    assessments: list[HypothesisAssessment] = []
    witnesses: list[Witness] = []
    coverage: dict[str, set[str]] = {assay.id: set() for assay in manifest.candidates}
    for hypothesis in manifest.hypotheses:
        differences = tuple(a.id for a in manifest.assays if
                            signatures[hypothesis.id, a.id] != signatures[expected.id, a.id])
        assessments.append(HypothesisAssessment(
            hypothesis_id=hypothesis.id, alleles=hypothesis.alleles,
            equivalent_to_expected=not differences,
            distinguishing_existing_assays=differences,
        ))
        if hypothesis.id == expected.id or differences:
            continue
        resolving = tuple(a.id for a in manifest.candidates if
                          signatures[hypothesis.id, a.id] != signatures[expected.id, a.id])
        for assay_id in resolving:
            coverage[assay_id].add(hypothesis.id)
        missing = sorted({
            allele_id for allele_id in hypothesis.alleles
            if any(by_pair[allele_id, a.id].signal_id is None for a in manifest.assays)
        })
        reason = (
            "Both hypotheses yield the same sequence-presence observations in every existing assay. "
            "At least one alternative allele emits no modeled signal in an existing assay: "
            + ", ".join(missing) + ". The surviving signal does not establish the missing allele's state."
            if missing else
            "Both hypotheses yield the same sequence-presence observations in every existing assay. "
            "The chosen readout does not distinguish the declared alleles; matching readouts do not "
            "establish matching full genomic states."
        )
        witnesses.append(Witness(
            hypothesis_id=hypothesis.id, expected_alleles=expected.alleles,
            alternative_alleles=hypothesis.alleles, explanation=reason,
            resolving_candidate_assays=resolving,
        ))
    notices = [Notice(
        code="FINITE_HYPOTHESIS_SPACE",
        message="Absence of a counterexample among declared hypotheses is not evidence of biological completeness.",
    ), Notice(
        code="NO_EXPERIMENTAL_VALIDATION",
        message="This alpha has software tests, not measured assay sensitivity or empirical validation.",
    )]
    if manifest.reference.synthetic:
        notices.append(Notice(code="SYNTHETIC_REFERENCE", message="This reference is synthetic demonstration data."))
    for assay in assays:
        left = manifest.reference.sequence[assay.left_primer.start:assay.left_primer.end]
        right = manifest.reference.sequence[assay.right_primer.start:assay.right_primer.end]
        if len(find_all(manifest.reference.sequence, left)) > 1 or len(
            find_all(manifest.reference.sequence, right)
        ) > 1:
            notices.append(Notice(
                code="LOCAL_PRIMER_MULTIMATCH", related_ids=(assay.id,),
                message="An annotated primer-site sequence has multiple exact local matches. "
                        "Only the annotated sites are modeled; genome-wide specificity is not assessed.",
            ))
        if not signatures[expected.id, assay.id]:
            notices.append(Notice(
                code="NO_EXPECTED_SIGNAL", related_ids=(assay.id,),
                message="The expected hypothesis yields no signal in this assay under the declared model. "
                        "A real negative result is also compatible with experimental failure.",
            ))
    if all(not signatures[expected.id, a.id] for a in manifest.assays):
        notices.append(Notice(
            code="NO_BASELINE_POSITIVE_SIGNAL",
            message="None of the existing assays yields expected signal; this design cannot positively confirm the expectation.",
        ))
    ref = manifest.reference
    result = Analysis(
        package_version=__version__, model_version=MODEL_VERSION,
        manifest_sha256=digest(manifest.model_dump(mode="json")),
        conclusion="ambiguity_demonstrated" if witnesses else "distinguishable_only_within_declared_model",
        expected_hypothesis=expected.id, expected_alleles=expected.alleles,
        reference=ReferenceSummary(
            name=ref.name, length=len(ref.sequence), sequence_sha256=hashlib.sha256(ref.sequence.encode()).hexdigest(),
            synthetic=ref.synthetic, assembly=ref.assembly, contig=ref.contig, genomic_start=ref.genomic_start,
        ),
        assays=manifest.assays, candidates=manifest.candidates, assumptions=ASSUMPTIONS,
        notices=tuple(notices), allele_observations=observations,
        hypothesis_observations=tuple(h_observations), hypotheses=tuple(assessments),
        witnesses=tuple(witnesses),
        plan=plan_panel(manifest.candidates, coverage, {w.hypothesis_id for w in witnesses}),
    )
    return seal(result)
