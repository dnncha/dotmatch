"""Pure, deterministic analysis; the report layer never changes scientific conclusions."""
from __future__ import annotations

import hashlib

from ._version import EXACT_MODEL_VERSION, __version__
from .io import InputError, digest, seal
from .models import (
    AlleleEvidence, Analysis, HypothesisAssessment, HypothesisObservation, Manifest, Notice,
    ReferenceSummary, Witness,
)
from .observations import observe_allele
from .planner import plan_panel
from .sequence import apply_edits, find_all
from .exact import EvidenceBudget, observe_exact

MAX_HYPOTHESIS_SIGNAL_REFERENCES = 100_000

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
    # Frozen models can still be constructed/copied without validation by callers.
    # Revalidate at the public computation boundary, including nested instances.
    manifest = Manifest.model_validate(manifest)
    assays = manifest.assays + manifest.candidates
    sequences = {a.id: apply_edits(manifest.reference.sequence, a.edits) for a in manifest.alleles}
    sequence_digests = {key: hashlib.sha256(seq.encode()).hexdigest() for key, seq in sequences.items()}
    budget = EvidenceBudget()
    expected = next(h for h in manifest.hypotheses if h.id == manifest.expected_hypothesis)
    observations = tuple(
        (observe_exact(manifest.reference.sequence, allele, assay,
                       edited_sequence=sequences[allele.id], budget=budget)
         if manifest.observation_model == EXACT_MODEL_VERSION else
         observe_allele(manifest.reference.sequence, allele, assay))
        for allele in manifest.alleles for assay in assays
    )
    by_pair = {(o.allele_id, o.assay_id): o for o in observations}
    h_observations: list[HypothesisObservation] = []
    signatures: dict[tuple[str, str], tuple[str, ...]] = {}
    signal_references = 0
    for hypothesis in manifest.hypotheses:
        for assay in assays:
            ids = tuple(sorted({
                signal for allele_id in hypothesis.alleles
                for signal in by_pair[allele_id, assay.id].signal_ids
            }))
            signal_references += len(ids)
            if signal_references > MAX_HYPOTHESIS_SIGNAL_REFERENCES:
                raise InputError(
                    "hypothesis signal evidence exceeds 100,000 references; narrow the hypothesis "
                    "space, use more specific primers or split the manifest. No partial result was returned."
                )
            signatures[hypothesis.id, assay.id] = ids
            h_observations.append(HypothesisObservation(
                hypothesis_id=hypothesis.id, assay_id=assay.id, signal_ids=ids
            ))
    assessments: list[HypothesisAssessment] = []
    witnesses: list[Witness] = []
    coverage: dict[str, set[str]] = {assay.id: set() for assay in manifest.candidates}
    expected_state = tuple(sorted(sequences[a] for a in expected.alleles))
    for hypothesis in manifest.hypotheses:
        same_state = tuple(sorted(sequences[a] for a in hypothesis.alleles)) == expected_state
        differences = tuple(a.id for a in manifest.assays if
                            signatures[hypothesis.id, a.id] != signatures[expected.id, a.id])
        assessments.append(HypothesisAssessment(
            hypothesis_id=hypothesis.id, alleles=hypothesis.alleles,
            equivalent_to_expected=not differences,
            distinguishing_existing_assays=differences,
            same_local_genomic_state_as_expected=same_state, description=hypothesis.description,
        ))
        if same_state or differences:
            continue
        resolving = tuple(a.id for a in manifest.candidates if
                          signatures[hypothesis.id, a.id] != signatures[expected.id, a.id])
        for assay_id in resolving:
            coverage[assay_id].add(hypothesis.id)
        missing = sorted({
            allele_id for allele_id in hypothesis.alleles
            if any(not by_pair[allele_id, a.id].signal_ids for a in manifest.assays)
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
    if manifest.observation_model == EXACT_MODEL_VERSION:
        assumptions = (ASSUMPTIONS[0],
            "Exact matching uses all inward-facing heteroprimer products in both orientations of each "
            "reconstructed local allele. New and rescued exact sites are included. Single-primer "
            "products, mismatches, off-window sites and thermodynamic efficiencies are not modeled.",
            *ASSUMPTIONS[2:])
        for observation in observations:
            if len(observation.products) > 1:
                notices.append(Notice(code="MULTIPLE_LOCAL_PRODUCTS", related_ids=(observation.assay_id, observation.allele_id),
                    message="All eligible exact local products contribute to the signal set; none was silently chosen or dropped."))
    else:
        assumptions = ASSUMPTIONS
        notices.append(Notice(code="LEGACY_REPRESENTATION_DEPENDENCE", message=
            "Original-site eligibility depends on edit notation. Use the explicitly selected exact-local v2 "
            "model for sequence-invariant rematching; neither model predicts real PCR efficiency."))
    identical = tuple(h.hypothesis_id for h in assessments
                      if h.hypothesis_id != expected.id and h.same_local_genomic_state_as_expected)
    if identical:
        notices.append(Notice(code="IDENTICAL_LOCAL_STATES_EXCLUDED", related_ids=identical,
            message="These hypotheses encode the same local diploid sequence state as the expectation. "
                    "Renaming or re-encoding a state does not create a biological counterexample."))
    if manifest.generation is not None:
        g = manifest.generation
        notices.append(Notice(code="GENERATED_DELETION_HYPOTHESES", message=
            f"Generation metadata declares {g.valid_deletions} grid deletions, {g.added_hypotheses} added "
            f"hypotheses and {g.deduplicated_states} duplicate states. Counts are not biological frequencies."))
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
                message=("An annotated primer sequence has multiple local matches. "
                         + ("All eligible exact heteroprimer products are modeled. "
                            if manifest.observation_model == EXACT_MODEL_VERSION else
                            "Only annotated original sites are modeled. ")
                         + "Genome-wide specificity is not assessed."),
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
        package_version=__version__, model_version=manifest.observation_model,
        manifest_sha256=digest(manifest.model_dump(mode="json")),
        conclusion="ambiguity_demonstrated" if witnesses else "distinguishable_only_within_declared_model",
        expected_hypothesis=expected.id, expected_alleles=expected.alleles,
        reference=ReferenceSummary(
            name=ref.name, length=len(ref.sequence), sequence_sha256=hashlib.sha256(ref.sequence.encode()).hexdigest(),
            synthetic=ref.synthetic, assembly=ref.assembly, contig=ref.contig, genomic_start=ref.genomic_start,
        ),
        assays=manifest.assays, candidates=manifest.candidates, assumptions=assumptions,
        notices=tuple(notices),
        generation=manifest.generation,
        allele_evidence=tuple(AlleleEvidence(allele_id=a.id, description=a.description,
                             edits=a.edits, sequence_length=len(sequences[a.id]),
                             sequence_sha256=sequence_digests[a.id]) for a in manifest.alleles),
        allele_observations=observations,
        hypothesis_observations=tuple(h_observations), hypotheses=tuple(assessments),
        witnesses=tuple(witnesses),
        plan=plan_panel(manifest.candidates, coverage, {w.hypothesis_id for w in witnesses}),
    )
    return seal(result)
