import sys, time, random, json, statistics, platform
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
from dotmatch.core import Matcher, MATCH_UNIQUE
from dotmatch.sensitivity import compare_hamming_policies

rng = random.Random(63104)
targets = ["".join(rng.choices("ACGT", k=20)) for _ in range(4000)]
reads = []
for i in range(100000):
    s = targets[rng.randrange(len(targets))]
    if i % 4 == 0:
        j = rng.randrange(20)
        s = s[:j] + rng.choice("ACGT") + s[j + 1 :]
    reads.append(s)
batches = [reads[i : i + 4096] for i in range(0, len(reads), 4096)]


def sequential(m):
    out = []
    for chunk in batches:
        modes = [
            m.assign_exact(chunk),
            m.assign_hamming(chunk, k=1, policy="radius"),
            m.assign_hamming(chunk, k=1, policy="best"),
        ]
        out.extend(
            tuple(
                (r.status, r.target_index if r.status == MATCH_UNIQUE else -1)
                for r in row
            )
            for row in zip(*modes)
        )
    return out


def fused(m):
    out = []
    for chunk in batches:
        out.extend(
            tuple(
                (getattr(row, mode).status, getattr(row, mode).target_index)
                for mode in ("exact", "radius_k1", "best_k1")
            )
            for row in compare_hamming_policies(m, chunk)
        )
    return out


measurements = {name: [] for name in ["sequential_three_queries", "fused_one_query"]}
with Matcher(targets) as m:
    sequential(m)
    fused(m)
    for i in range(5):
        methods = [("sequential_three_queries", sequential), ("fused_one_query", fused)]
        if i % 2:
            methods.reverse()
        results = []
        for name, method in methods:
            start = time.perf_counter()
            result = method(m)
            measurements[name].append(time.perf_counter() - start)
            results.append(result)
        assert results[0] == results[1]
result = {
    "kind": "local_synthetic_engineering_measurement",
    "seed": 63104,
    "targets": 4000,
    "reads": 100000,
    "length": 20,
    "batch_size": 4096,
    "repeats": 5,
    "python": sys.version.split()[0],
    "platform": platform.platform(),
    "seconds": measurements,
    "median_seconds": {k: statistics.median(v) for k, v in measurements.items()},
    "all_policy_outcomes_and_unique_ids_equal": True,
    "limits": "One process, prebuilt warmed native index, synthetic reads, allocation and Python projection included; excludes file I/O and index construction. Not a competitor or end-to-end FASTQ benchmark.",
}
result["median_ratio"] = (
    result["median_seconds"]["sequential_three_queries"]
    / result["median_seconds"]["fused_one_query"]
)
(ROOT / "benchmarks/raw/policy_sensitivity_local.json").write_text(
    json.dumps(result, indent=2) + "\n"
)
print(json.dumps(result, indent=2))
