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
