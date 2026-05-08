import dotmatch


def test_dotmatch_exports_public_api():
    assert dotmatch.distance("ACGT", "AGGT") == 1
    assert dotmatch.distance_leq("ACGT", "AGGT", 1)
    assert hasattr(dotmatch, "assign")


def test_dotmatch_exports_literal_byte_alphabet_policy():
    assert dotmatch.alphabet_policy() == (
        "literal-byte; A/C/G/T/N/IUPAC symbols are ordinary byte symbols; no wildcard expansion"
    )
    assert dotmatch.distance("N", "A") == 1
    assert dotmatch.distance("R", "A") == 1
    assert dotmatch.distance("R", "R") == 0
