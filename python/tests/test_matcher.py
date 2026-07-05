import dotmatch
import quickdna
import pytest


def test_matcher_assign_matches_scan_api():
    reads = ["ACGT", "ACGC", "TTTT"]
    targets = ["ACGT", "AGGT", "ACGA"]

    matcher = dotmatch.Matcher(targets)
    indexed = matcher.assign(reads, k=1)
    scan = dotmatch.assign(reads, targets, k=1)

    assert indexed == scan


def test_matcher_assign_with_stats_exposes_candidate_counts():
    matcher = dotmatch.Matcher(["AAAAAAAA", "CCCCCCCC", "GGGGGGGG", "TTTTTTTT"])

    results, stats = matcher.assign_with_stats(["AAAAAAAT", "CCCCCCCA"], k=1)

    assert [r.status for r in results] == [dotmatch.MATCH_UNIQUE, dotmatch.MATCH_UNIQUE]
    assert stats.candidates_considered == stats.candidates_verified
    assert 0 < stats.candidates_verified < 8


def test_matcher_exact_uses_direct_lookup_and_supports_ascii_fold():
    matcher = dotmatch.Matcher(["ACGT", "ACGT", "TTTT", "NNNN"])

    exact, stats = matcher.assign_exact_with_stats(["ACGT", "acgt", "NNNN", "GGGG"])
    folded, folded_stats = matcher.assign_exact_with_stats(["ACGT", "acgt", "NNNN", "nnnn"], ascii_fold=True)

    assert [r.status for r in exact] == [
        dotmatch.MATCH_AMBIGUOUS,
        dotmatch.MATCH_NONE,
        dotmatch.MATCH_UNIQUE,
        dotmatch.MATCH_NONE,
    ]
    assert folded[0].status == dotmatch.MATCH_AMBIGUOUS
    assert folded[1].status == dotmatch.MATCH_AMBIGUOUS
    assert folded[2].target_index == 3
    assert folded[3].target_index == 3
    assert stats.candidates_verified == 4
    assert folded_stats.candidates_verified == 6


def test_matcher_hamming_matches_levenshtein_for_equal_length_substitutions():
    reads = ["ACGTACGT", "ACGTACGA", "TTTTTTTT", "ACGTACGG"]
    targets = ["ACGTACGT", "ACGTACGA", "CCCCCCCC", "GGGGGGGG"]

    matcher = dotmatch.Matcher(targets)
    hamming, stats = matcher.assign_hamming_with_stats(reads, k=1, policy="best")
    levenshtein = matcher.assign(reads, k=1, policy="best")

    assert hamming == levenshtein
    assert 0 < stats.candidates_verified < len(reads) * len(targets)


def test_top_level_hamming_and_exact_helpers():
    targets = ["ACGT", "TTTT"]

    assert dotmatch.assign_hamming(["ACGA"], targets, k=1, policy="best")[0].target_index == 0
    assert dotmatch.assign_exact(["acgt"], targets, ascii_fold=True)[0].target_index == 0


def test_assign_dataframe_metric_hamming_routes_to_fixed_length_kernel():
    pd = pytest.importorskip("pandas")
    reads = pd.Series(["ACGT", "ACGA", "ACGTT"], index=["exact", "sub", "indel"])
    targets = pd.DataFrame({"id": ["a"], "seq": ["ACGT"]})

    hamming = dotmatch.assign_dataframe(reads, targets, k=1, metric="hamming", policy="best")
    levenshtein = dotmatch.assign_dataframe(reads, targets, k=1, metric="levenshtein", policy="best")

    assert hamming.loc[0, "status_name"] == "unique"
    assert hamming.loc[1, "status_name"] == "unique"
    assert hamming.loc[2, "status_name"] == "none"
    assert levenshtein.loc[2, "status_name"] == "unique"


def test_quickdna_compatibility_exports_matcher():
    assert quickdna.Matcher(["ACGT"]).assign(["ACGT"], k=0)[0].status == quickdna.MATCH_UNIQUE


def test_assign_defaults_to_radius_safe_ambiguity_policy():
    result = dotmatch.assign(["ACGT"], ["ACGT", "AGGT", "ACGA"], k=1)[0]

    assert result.status == dotmatch.MATCH_AMBIGUOUS
    assert result.target_index == 0
    assert result.best_distance == 0
    assert result.second_best_distance == 1
    assert result.match_count == 3


def test_assign_can_opt_into_best_distance_policy():
    result = dotmatch.assign(["ACGT"], ["ACGT", "AGGT", "ACGA"], k=1, policy="best")[0]

    assert result.status == dotmatch.MATCH_UNIQUE
    assert result.target_index == 0
    assert result.best_distance == 0
    assert result.match_count == 3


def test_assign_posterior_accepts_confident_exact_call():
    result = dotmatch.assign_posterior("ACGT", ["ACGT", "ACGA"], "IIII", min_posterior=0.95)

    assert result.status == dotmatch.MATCH_UNIQUE
    assert result.target_index == 0
    assert result.posterior > 0.99
    assert result.second_posterior < 0.01


def test_assign_posterior_rejects_low_quality_tie():
    result = dotmatch.assign_posterior("ACGT", ["ACGT", "ACGA"], "III#", min_posterior=0.95)

    assert result.status == dotmatch.MATCH_AMBIGUOUS
    assert result.target_index == 0
    assert result.posterior < 0.95


def test_assign_posterior_uses_quality_to_make_confident_call():
    result = dotmatch.assign_posterior("ACGT", ["ACGT", "AGGT"], "IIII")

    assert result.status == dotmatch.MATCH_UNIQUE
    assert result.target_index == 0
    assert result.posterior > 0.99
    assert result.second_posterior < 0.01


def test_assign_posterior_reports_ambiguity_when_quality_is_weak():
    result = dotmatch.assign_posterior("ACGT", ["ACGT", "AGGT"], "####")

    assert result.status == dotmatch.MATCH_AMBIGUOUS
    assert result.target_index == 0
    assert result.posterior < 0.95
    assert len(result.posteriors) == 2


def test_assign_posterior_accepts_non_uniform_priors():
    result = dotmatch.assign_posterior("ACGT", ["ACGT", "AGGT"], "####", priors=[0.01, 0.99])

    assert result.status == dotmatch.MATCH_UNIQUE
    assert result.target_index == 1


def test_assign_posterior_marks_length_mismatch_invalid():
    result = dotmatch.assign_posterior("ACGT", ["ACGT", "ACGTA"], "IIII")

    assert result.status == dotmatch.MATCH_INVALID
    assert result.target_index == -1
    assert result.posteriors == ()


def test_assign_posterior_validates_quality_length():
    try:
        dotmatch.assign_posterior("ACGT", ["ACGT"], "III")
    except ValueError as exc:
        assert "same length" in str(exc)
    else:
        raise AssertionError("expected quality length validation")
