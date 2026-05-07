import dotmatch


def test_dotmatch_exports_public_api():
    assert dotmatch.distance("ACGT", "AGGT") == 1
    assert dotmatch.distance_leq("ACGT", "AGGT", 1)
    assert hasattr(dotmatch, "assign")
