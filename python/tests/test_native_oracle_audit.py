"""Seeded scalar dynamic-programming oracle independent of the native kernels."""
import random
import pytest
import dotmatch


def levenshtein(a, b):
    previous = list(range(len(b) + 1))
    for i, x in enumerate(a, 1):
        current = [i]
        for j, y in enumerate(b, 1):
            current.append(min(current[-1] + 1, previous[j] + 1, previous[j-1] + (x != y)))
        previous = current
    return previous[-1]


@pytest.mark.parametrize('length', [0, 1, 8, 20, 31, 32, 33, 63, 64, 65, 127, 128, 129])
def test_distance_boundaries_against_scalar_dp(length):
    rng = random.Random(8849 + length)
    for _ in range(8):
        a = ''.join(rng.choices('ACGTNacgt', k=length))
        b = list(a)
        if b:
            b[rng.randrange(len(b))] = rng.choice('ACGTN')
        if rng.randrange(2): b.append('A')
        b = ''.join(b)
        expected = levenshtein(a, b)
        assert dotmatch.distance(a, b) == expected
        for k in [0, 1, 2, 3, 128]:
            assert dotmatch.distance_leq(a, b, k) == (expected <= k)


@pytest.mark.parametrize('length', [8, 20, 32, 65])
def test_indexed_and_unindexed_policies_agree_with_scalar_oracle(length):
    rng = random.Random(length)
    first = ''.join(rng.choices('ACGT', k=length))
    targets = [first, first]  # duplicated sequence is intentionally ambiguous
    for i in range(8):
        value = list(first)
        for _ in range(i % 4): value[rng.randrange(length)] = rng.choice('ACGT')
        targets.append(''.join(value))
    reads = targets + [first[:-1], 'N' * length, first + 'A']
    with dotmatch.Matcher(targets) as matcher:
        for k in range(3):
            for policy in ['radius', 'best']:
                outputs = [dotmatch.assign(reads, targets, k, policy), matcher.assign(reads, k, policy)]
                for output in outputs:
                    for read, result in zip(reads, output):
                        distances = [levenshtein(read, target) for target in targets]
                        nearest = min(distances)
                        candidates = [i for i, d in enumerate(distances) if d <= k and (policy == 'radius' or d == nearest)]
                        assert result.status == (0 if not candidates else 1 if len(candidates) == 1 else 2)
                        if len(candidates) == 1: assert result.target_index == candidates[0]
