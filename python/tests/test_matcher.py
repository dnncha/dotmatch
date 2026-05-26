import dotmatch
import quickdna


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
