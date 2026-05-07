import quickdna


def test_distance_and_threshold():
    assert quickdna.distance("ACGT", "AGGT") == 1
    assert quickdna.distance_leq("ACGT", "AGGT", 1)
    assert not quickdna.distance_leq("ACGT", "AGGT", 0)


def test_assign_unique_ambiguous_none():
    results = quickdna.assign(["ACGT", "ACGC", "TTTT"], ["ACGT", "AGGT", "ACGA"], k=1)

    assert results[0].status == quickdna.MATCH_UNIQUE
    assert results[0].target_index == 0
    assert results[0].best_distance == 0
    assert results[0].match_count == 3

    assert results[1].status == quickdna.MATCH_AMBIGUOUS
    assert results[1].best_distance == 1
    assert results[1].match_count == 2

    assert results[2].status == quickdna.MATCH_NONE
    assert results[2].target_index == -1


def test_assign_rejects_negative_k():
    try:
        quickdna.assign(["ACGT"], ["ACGT"], k=-1)
    except ValueError:
        pass
    else:
        raise AssertionError("negative k should raise ValueError")
