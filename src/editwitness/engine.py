"""Pure, deterministic analysis; the report layer never changes scientific conclusions."""
from __future__ import annotations

import hashlib
from typing import Literal

from ._version import EXACT_MODEL_VERSION, __version__
from .exact import ProductBudget, observe_exact, signal_ids
from .io import digest, seal
from .models import (
    Analysis, HypothesisAssessment, HypothesisObservation, Manifest, Notice,
    ReferenceSummary, Witness, validated_manifest,
)
from .observations import observe_allele
from .planner import plan_panel
from .sequence import apply_edits, find_all

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
    manifest = validated_manifest(manifest)
    assays = manifest.assays + manifest.candidates
    reconstructed = {a.id: apply_edits(manifest.reference.sequence, a.edits) for a in manifest.alleles}
    exact = manifest.observation_model == EXACT_MODEL_VERSION
    budget = ProductBudget()
    expected = next(h for h in manifest.hypotheses if h.id == manifest.expected_hypothesis)
    observations = tuple(
        (observe_exact(manifest.reference.sequence, allele, assay,
                       edited_sequence=reconstructed[allele.id], budget=budget)
         if exact else observe_allele(manifest.reference.sequence, allele, assay))
        for allele in manifest.alleles for assay in assays
    )
    by_pair = {(o.allele_id, o.assay_id): o for o in observations}
    h_observations: list[HypothesisObservation] = []
    signatures: dict[tuple[str, str], tuple[str, ...]] = {}
    for hypothesis in manifest.hypotheses:
        for assay in assays:
            ids = tuple(sorted({
                signal for allele_id in hypothesis.alleles
                for signal in signal_ids(by_pair[allele_id, assay.id])
            }))
            signatures[hypothesis.id, assay.id] = ids
            h_observations.append(HypothesisObservation(
                hypothesis_id=hypothesis.id, assay_id=assay.id, signal_ids=ids
            ))
    expected_genotype = tuple(sorted(reconstructed[a] for a in expected.alleles))
    genotypes = {h.id: tuple(sorted(reconstructed[a] for a in h.alleles)) for h in manifest.hypotheses}
    representatives: dict[tuple[str, ...], str] = {}
    for h in sorted(manifest.hypotheses, key=lambda h: h.id):
        representatives.setdefault(genotypes[h.id], h.id)
    representatives[expected_genotype] = expected.id
    same_genotype = {identifier: genotype == expected_genotype for identifier, genotype in genotypes.items()}
    distinct_alternatives = len(representatives)-1
    assessments: list[HypothesisAssessment] = []
    witnesses: list[Witness] = []
    coverage: dict[str, set[str]] = {assay.id: set() for assay in manifest.candidates}
    for hypothesis in manifest.hypotheses:
        differences = tuple(a.id for a in manifest.assays if
                            signatures[hypothesis.id, a.id] != signatures[expected.id, a.id])
        assessments.append(HypothesisAssessment(
            hypothesis_id=hypothesis.id, alleles=hypothesis.alleles,
            equivalent_to_expected=not differences,
            same_local_genotype_as_expected=same_genotype[hypothesis.id],
            representative_hypothesis=representatives[genotypes[hypothesis.id]],
            distinguishing_existing_assays=differences,
        ))
        if same_genotype[hypothesis.id] or differences or representatives[genotypes[hypothesis.id]] != hypothesis.id:
            continue
        resolving = tuple(a.id for a in manifest.candidates if
                          signatures[hypothesis.id, a.id] != signatures[expected.id, a.id])
        for assay_id in resolving:
            coverage[assay_id].add(hypothesis.id)
        missing = sorted({
            allele_id for allele_id in hypothesis.alleles
            if any(not signal_ids(by_pair[allele_id, a.id]) for a in manifest.assays)
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
    aliases = tuple(h.id for h in manifest.hypotheses if h.id != expected.id and same_genotype[h.id])
    if aliases:
        notices.append(Notice(code="SAME_LOCAL_GENOTYPE", related_ids=aliases,
                              message="These hypotheses reconstruct the same unordered local diploid sequences "
                                      "as the expectation. They are not distinct counterexamples."))
    duplicate_ids = tuple(h.id for h in manifest.hypotheses
                          if not same_genotype[h.id] and representatives[genotypes[h.id]] != h.id)
    if duplicate_ids:
        notices.append(Notice(code="DUPLICATE_LOCAL_GENOTYPES", related_ids=duplicate_ids,
                              message="Sequence-identical alternative genotypes share one representative witness. "
                                      "Aliases do not inflate counterexample counts or panel coverage."))
    if not distinct_alternatives:
        notices.append(Notice(code="NO_DISTINCT_ALTERNATIVES",
                              message="No different local genotype was supplied; no assay discrimination was tested."))
    if exact:
        notices.append(Notice(code="EXACT_MATCH_MODEL",
                              message="Both inward heteroprimer orientations are rematched on edited local sequences. "
                                      "Exact matching is not PCR prediction; mismatches, same-primer products and "
                                      "sites outside this window are not modeled."))
        multiple = tuple(sorted({o.assay_id for o in observations if len(o.products) > 1}))
        if multiple:
            notices.append(Notice(code="MULTIPLE_LOCAL_PRODUCTS", related_ids=multiple,
                                  message="Some allele/assay combinations yield several modeled products. "
                                          "All sequence signals are retained; no dominant product is assumed."))
    else:
        notices.append(Notice(code="LEGACY_ORIGINAL_SITE_MODEL",
                              message="Original-site disruption is representation dependent. Use compare-models "
                                      "to inspect sensitivity to exact local sequence rematching."))
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
                        + ("All exact local heteroprimer products are modeled; genome-wide specificity is not assessed."
                         if exact else "Only the annotated sites are modeled; genome-wide specificity is not assessed."),
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
    baseline_informative = any(signatures[expected.id, a.id] for a in manifest.assays)
    conclusion: Literal["no_distinct_alternatives", "baseline_uninformative", "ambiguity_demonstrated", "distinguishable_only_within_declared_model"] = (
        "no_distinct_alternatives" if not distinct_alternatives else
        "baseline_uninformative" if not baseline_informative else
        "ambiguity_demonstrated" if witnesses else "distinguishable_only_within_declared_model"
    )
    assumptions = list(ASSUMPTIONS)
    if exact:
        assumptions[1] = (
            "Exact primer matches are searched on each reconstructed local allele, in both inward-facing "
            "heteroprimer orientations. Every in-bounds product contributes a sequence signal. "
            "Mismatches, same-primer amplification and nonlocal binding sites are not modeled."
        )
    ref = manifest.reference
    result = Analysis(
        package_version=__version__, model_version=manifest.observation_model,
        manifest_sha256=digest(manifest.model_dump(mode="json")),
        conclusion=conclusion, distinct_alternatives=distinct_alternatives,
        alleles=manifest.alleles, generation=manifest.generation,
        expected_hypothesis=expected.id, expected_alleles=expected.alleles,
        reference=ReferenceSummary(
            name=ref.name, length=len(ref.sequence), sequence_sha256=hashlib.sha256(ref.sequence.encode()).hexdigest(),
            synthetic=ref.synthetic, assembly=ref.assembly, contig=ref.contig, genomic_start=ref.genomic_start,
        ),
        assays=manifest.assays, candidates=manifest.candidates, assumptions=tuple(assumptions),
        notices=tuple(notices), allele_observations=observations,
        hypothesis_observations=tuple(h_observations), hypotheses=tuple(assessments),
        witnesses=tuple(witnesses),
        plan=plan_panel(manifest.candidates, coverage, {w.hypothesis_id for w in witnesses}),
    )
    return seal(result)
