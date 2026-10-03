"""Held-out essential-gene ranking on two pinned public HAP1 screens.

Run BAGEL2 with external, gene-stratified training folds. No evaluation label
enters its fold's training files. Fit is unsupervised. All methods see the same
complete-four-guide cohort and samples; preprocessing is method-specific.
"""

import argparse
from collections import Counter
import csv
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

import numpy as np


COMMIT = "53388adbb4fb0931e5c9dda135502be19e4555f0"
HASHES = {
    "BAGEL.py": "89f9ba4115513caf9838ce40ab4ca8d47e5e2bc82cd128f356185997ddf9a6c7",
    "CEGv2.txt": "1b1df46f829e2c47dc8a4c164f60f5f1dc67513bafba7d8b9609421e0c75d27e",
    "NEGv1.txt": "18a6b37a8f182a1bf501548e6b28fe4fcef02fcd296f873ef500a21eeb433739",
    "reads_hap1.txt": "7638bc6237cf8a3e2302fcd969d014b041819669b37f2f87bc1f7cfbc45ca8a7",
    "pipeline-script-example/HAP1-TKOv3-EXAMPLE.txt": "9c8e3c81efef6b1761a28ec54813f5c3c8d596466f8dc85caa49d673817681a1",
}
SCREENS = {
    "hap1-tko": ("reads_hap1.txt", ["HAP1_T0", "HAP1_T18A", "HAP1_T18B", "HAP1_T18C"]),
    "hap1-tkov3": ("pipeline-script-example/HAP1-TKOv3-EXAMPLE.txt", ["T0_-", "T18_A", "T18_B", "T18_C"]),
}
METHODS = ("fit_beta", "fit_wald_z", "fit_finite_depletion_p", "fit_control_beta",
           "fit_control_depletion_p", "bagel2", "mean_log2fc")


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def ranking_metrics(labels, scores):
    """Threshold-block AP and ROC AUC; score ties never use gene order."""
    labels = np.asarray(labels)
    scores = np.asarray(scores, dtype=float)
    if labels.ndim != 1 or scores.shape != labels.shape or not np.isfinite(scores).all():
        raise ValueError("Labels and finite scores must be matching 1D arrays")
    if not np.isin(labels, [0, 1]).all() or len(np.unique(labels)) != 2:
        raise ValueError("Both binary label classes are required")
    order = np.argsort(-scores, kind="stable")
    ranked = labels[order]
    ends = np.r_[np.flatnonzero(np.diff(scores[order]) != 0), len(scores) - 1]
    tp = np.cumsum(ranked)[ends]
    fp = ends + 1 - tp
    recall = tp / labels.sum()
    precision = tp / (ends + 1)
    fpr = np.r_[0., fp / (len(labels) - labels.sum())]
    tpr = np.r_[0., recall]
    return {
        "average_precision": float(np.sum(np.diff(np.r_[0., recall]) * precision)),
        "roc_auc": float(np.sum(np.diff(fpr) * (tpr[1:] + tpr[:-1]) / 2)),
        "recall_at_95pct_precision": float(np.max(recall[precision >= .95], initial=0)),
    }


def heldout_folds(essentials, nonessentials, *, seed=1729, folds=5):
    if folds < 2 or min(len(essentials), len(nonessentials)) < folds:
        raise ValueError("At least two folds and one gene from each class per fold are required")
    if set(essentials) & set(nonessentials):
        raise ValueError("Essential and nonessential labels overlap")
    rng = np.random.default_rng(seed)
    assignment = {}
    for genes in (essentials, nonessentials):
        shuffled = np.array(sorted(genes))
        rng.shuffle(shuffled)
        assignment.update({gene: i % folds for i, gene in enumerate(shuffled)})
    return assignment


def write_rows(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def summarize(rows, *, seed=1729, bootstraps=500):
    metrics = []
    for dataset in sorted({row["dataset"] for row in rows}):
        cohort = [row for row in rows if row["dataset"] == dataset]
        labels = np.array([int(row["essential"]) for row in cohort])
        positives, negatives = np.flatnonzero(labels == 1), np.flatnonzero(labels == 0)
        rng = np.random.default_rng(seed)
        draws = [np.r_[rng.choice(positives, len(positives)), rng.choice(negatives, len(negatives))]
                 for _ in range(bootstraps)]
        for method in METHODS:
            scores = np.array([float(row[method]) for row in cohort])
            result = ranking_metrics(labels, scores)
            bootstrap_ap = [ranking_metrics(labels[draw], scores[draw])["average_precision"] for draw in draws]
            calls = None
            if method in ("fit_finite_depletion_p", "fit_control_depletion_p"):
                lane = "fit_control" if method == "fit_control_depletion_p" else "fit"
                calls = np.array([float(row[lane + "_depletion_fdr"]) <= .05
                                  and float(row[lane + "_beta"]) > 0 for row in cohort])
            metrics.append({
                "dataset": dataset, "method": method, "genes": len(cohort),
                "essentials": len(positives), "nonessentials": len(negatives), **result,
                "ap_ci_low": float(np.quantile(bootstrap_ap, .025)),
                "ap_ci_high": float(np.quantile(bootstrap_ap, .975)),
                "essential_call_rate_at_nominal_fdr_005": float(calls[positives].mean()) if calls is not None else "",
                "nonessential_call_rate_at_nominal_fdr_005": float(calls[negatives].mean()) if calls is not None else "",
            })
    return metrics


def report(scores_path, output_dir, *, seed=1729, bootstraps=500):
    with scores_path.open() as stream:
        rows = list(csv.DictReader(stream))
    metrics = summarize(rows, seed=seed, bootstraps=bootstraps)
    write_rows(output_dir / "metrics.csv", metrics)
    text = ["# CRISPRWorks Fit: held-out gene essentiality", "",
            "This evaluates ranking against CEGv2/NEGv1 annotations, not biological ground truth or calibrated FDR.",
            "Both screens are HAP1 from the Hart Lab BAGEL repository; they are not independent cell lines or labs.",
            "Five external stratified gene folds exclude each evaluation gene from BAGEL2's training labels.",
            "All methods receive complete four-guide labels, selected before normalization; incomplete labels and large control bins are excluded.",
            "Fit uses a T0 baseline and three T18 replicates. TKOv3 T3 and starvation samples are excluded.", "",
            "| Screen | Method | AP (95% conditional bootstrap CI) | ROC AUC | Recall at 95% precision | Essential / nonessential call rate at nominal 5% FDR |",
            "| --- | --- | --- | --- | --- | --- |"]
    for row in metrics:
        call_rates = "—" if row['essential_call_rate_at_nominal_fdr_005'] == "" else (
            f"{row['essential_call_rate_at_nominal_fdr_005']:.4f} / {row['nonessential_call_rate_at_nominal_fdr_005']:.4f}")
        text.append(f"| {row['dataset']} | {row['method']} | {row['average_precision']:.4f} "
                    f"({row['ap_ci_low']:.4f}–{row['ap_ci_high']:.4f}) | {row['roc_auc']:.4f} | "
                    f"{row['recall_at_95pct_precision']:.4f} | {call_rates} |")
    text += ["", "AP is threshold-block average precision; all ties enter together. Larger scores mean greater depletion.",
             f"Confidence intervals use {bootstraps} class-stratified gene bootstrap samples with seed {seed}, conditional on the fitted scores and training folds.",
             "They do not capture variability across screens or refitting the training sets. Recall at 95% precision is a retrospective ranking diagnostic, not a deployable FDR threshold.",
             "Fit scores are negative beta, negative Wald z, or negative finite depletion p. BAGEL2 is a supervised Bayes-factor classifier; mean log2FC is a simple baseline.",
             "BAGEL2 uses its own sum-read normalization and pseudocount 5, no network boost or multi-target correction; Fit uses median normalization and updated guide efficiency.",
             "The Fit control lane instead normalizes on nonessential training genes and draws its permutation background exclusively from their guides. Each scored gene is excluded from those controls.",
             "Call rates use negative fitted effects with directional permutation FDR <=0.05. They measure sensitivity and false-positive frequency in the annotated classes; the nonessential call rate is not the false discovery rate among all genome-wide hits.",
             "Finite tails repair tie/zero-p behavior. They do not change the fitted beta or Wald z and do not establish exchangeability of pooled pseudo-genes.",
             "", "Chronos 2.3 hit calling and JACKS, copy-number effects, other cell lines, multi-condition contrasts, and independent validation remain untested here.",
             "These measurements cannot support an overall state-of-the-art claim.", ""]
    (output_dir / "HIT-QUALITY.md").write_text("\n".join(text))
    return metrics


def run_command(command, log, provenance):
    started = time.perf_counter()
    with log.open("w") as stream:
        completed = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT, timeout=1800)
    provenance["runs"].append({"command": command, "wall_seconds": time.perf_counter() - started,
                               "exit_code": completed.returncode, "log": str(log)})
    if completed.returncode:
        raise RuntimeError(f"Benchmark command failed; see {log}")


def benchmark(bagel_root, output_dir, *, seed=1729, rounds=2):
    for name, expected in HASHES.items():
        if sha256(bagel_root / name) != expected:
            raise ValueError(f"Pinned source hash differs: {name}")
    from importlib.metadata import version
    from crisprworks_fit.kernels import resolved_engine
    provenance = {
        "schema_version": 1, "source_commit": COMMIT, "source_hashes": HASHES,
        "seed": seed, "folds": 5, "permutation_rounds": rounds, "runs": [],
        "python": platform.python_version(), "platform": platform.platform(),
        "versions": {name: version(name) for name in ("crisprworks-fit", "mageck2", "numpy", "scipy", "pandas", "scikit-learn", "joblib", "click")},
        "fit_engine": resolved_engine("native"),
        "source_urls": {name: f"https://raw.githubusercontent.com/hart-lab/bagel/{COMMIT}/{name}" for name in HASHES},
        "implementation_hashes": {str(path.relative_to(Path(__file__).parents[1])): sha256(path) for path in (
            Path(__file__), *(Path(__file__).parents[1] / "src" / "crisprworks_fit" / name for name in ("inference.py", "backend.py", "cli.py", "kernels.py", "_native.c")))},
    }
    labels = {}
    duplicates = []
    for name, value in (("CEGv2.txt", 1), ("NEGv1.txt", 0)):
        with (bagel_root / name).open() as stream:
            next(stream)
            for line in stream:
                gene = line.split("\t")[0].strip()
                if gene in labels:
                    if labels[gene] != value:
                        raise ValueError(f"Overlapping essential/nonessential label: {gene}")
                    duplicates.append(gene)
                labels[gene] = value
    provenance["duplicate_reference_rows_collapsed"] = duplicates
    all_scores = []
    provenance["screens"] = {}
    for dataset, (filename, samples) in SCREENS.items():
        directory = output_dir / dataset
        directory.mkdir(parents=True, exist_ok=True)
        with (bagel_root / filename).open() as stream:
            reader = csv.reader(stream, delimiter="\t")
            header = next(reader)
            source_rows = list(reader)
        guide_counts = Counter(row[1] for row in source_rows)
        selected = [row for row in source_rows if guide_counts[row[1]] == 4]
        columns = [header.index(sample) for sample in samples]
        counts = directory / "counts.tsv"
        with counts.open("w") as stream:
            stream.write("sgRNA\tGene\t" + "\t".join(samples) + "\n")
            for row in selected:
                stream.write("\t".join(row[:2] + [row[i] for i in columns]) + "\n")
        design = directory / "design.tsv"
        design.write_text("Samples\tbaseline\ttreatment\n" + "".join(
            f"{sample}\t1\t{int(i > 0)}\n" for i, sample in enumerate(samples)))
        prefix = directory / "fit"
        run_command([sys.executable, "-m", "crisprworks_fit", "mle", "-k", str(counts),
                     "-d", str(design), "-n", str(prefix), "--kernel", "native", "--seed", str(seed),
                     "--threads", "1", "--blas-threads", "1", "--permutation-round", str(rounds),
                     "--update-efficiency", "--genes-varmodeling", "1000", "--write-fit-details",
                     "--permutation-pvalues", "finite"], directory / "fit.stdout.log", provenance)
        details = json.loads(Path(str(prefix) + ".fit-details.json").read_text())["genes"]
        fit_manifest = json.loads(Path(str(prefix) + ".crisprworks.json").read_text())
        provenance["screens"][dataset] = {
            "counts_sha256": sha256(counts), "design_sha256": sha256(design),
            "guides": len(selected), "gene_labels": len(details), "samples": samples,
            "selection": "exactly four guides per gene before normalization; T0 plus three T18 samples",
            "missing_reference_genes": sorted(set(labels) - set(details)),
            "fit_manifest": fit_manifest,
        }
        essential = [gene for gene in details if labels.get(gene) == 1]
        nonessential = [gene for gene in details if labels.get(gene) == 0]
        folds = heldout_folds(essential, nonessential, seed=seed)
        fc_prefix = directory / "bagel"
        run_command([sys.executable, str(bagel_root / "BAGEL.py"), "fc", "-i", str(counts),
                     "-o", str(fc_prefix), "-c", samples[0]], directory / "bagel-fc.log", provenance)
        foldchange = Path(str(fc_prefix) + ".foldchange")
        mean_lfc = {}
        with foldchange.open() as stream:
            for row in csv.DictReader(stream, delimiter="\t"):
                mean_lfc.setdefault(row["GENE"], []).append(np.mean([float(row[sample]) for sample in samples[1:]]))
        bagel_scores = {}
        control_results = {}
        for fold in range(5):
            train_paths = []
            for name, genes in (("essential", essential), ("nonessential", nonessential)):
                train = sorted(gene for gene in genes if folds[gene] != fold)
                evaluation = {gene for gene in folds if folds[gene] == fold}
                if set(train) & evaluation:
                    raise RuntimeError("Evaluation genes leaked into training")
                path = directory / f"fold-{fold}-{name}.txt"
                path.write_text("GENE\n" + "\n".join(train) + "\n")
                train_paths.append(path)
            bf = directory / f"fold-{fold}.bf.tsv"
            run_command([sys.executable, str(bagel_root / "BAGEL.py"), "bf", "-i", str(foldchange),
                         "-o", str(bf), "-e", str(train_paths[0]), "-n", str(train_paths[1]),
                         "-c", ",".join(samples[1:]), "-NS"], directory / f"bagel-fold-{fold}.log", provenance)
            with bf.open() as stream:
                for row in csv.DictReader(stream, delimiter="\t"):
                    if folds.get(row["GENE"]) == fold:
                        bagel_scores[row["GENE"]] = float(row["BF"])
            # MAGeCK's control-gene input has no header; keep it separate from
            # BAGEL's header-bearing label file. Use only this fold's training
            # nonessentials for BOTH normalization and the pseudo-gene null.
            control_genes = directory / f"fold-{fold}-fit-controls.txt"
            control_genes.write_text("\n".join(sorted(gene for gene in nonessential if folds[gene] != fold)) + "\n")
            control_prefix = directory / f"fit-control-{fold}"
            run_command([sys.executable, "-m", "crisprworks_fit", "mle", "-k", str(counts),
                         "-d", str(design), "-n", str(control_prefix), "--kernel", "native", "--seed", str(seed),
                         "--threads", "1", "--blas-threads", "1", "--permutation-round", str(rounds),
                         "--update-efficiency", "--genes-varmodeling", "1000", "--write-fit-details",
                         "--permutation-pvalues", "finite", "--norm-method", "control", "--control-gene", str(control_genes)],
                        directory / f"fit-control-{fold}.stdout.log", provenance)
            control_details = json.loads(Path(str(control_prefix) + ".fit-details.json").read_text())["genes"]
            control_manifest = Path(str(control_prefix) + ".crisprworks.json")
            provenance["screens"][dataset].setdefault("control_fit_manifests", []).append(json.loads(control_manifest.read_text()))
            for gene in folds:
                if folds[gene] == fold:
                    control_results[gene] = control_details[gene]
        if set(bagel_scores) != set(folds):
            raise RuntimeError("BAGEL2 did not score every held-out reference gene")
        for gene in sorted(folds):
            result = details[gene]
            control = control_results[gene]
            all_scores.append({
                "dataset": dataset, "gene": gene, "essential": labels[gene], "fold": folds[gene],
                "fit_beta": -result["beta_estimate"][4], "fit_wald_z": -result["beta_zscore"][0],
                "fit_finite_depletion_p": -result["beta_permute_pval_neg"][0],
                "fit_depletion_fdr": result["beta_permute_pval_neg_fdr"][0],
                "fit_control_beta": -control["beta_estimate"][4],
                "fit_control_depletion_p": -control["beta_permute_pval_neg"][0],
                "fit_control_depletion_fdr": control["beta_permute_pval_neg_fdr"][0],
                "bagel2": bagel_scores[gene], "mean_log2fc": -float(np.mean(mean_lfc[gene])),
            })
        # Persist progress so a failure in a later screen retains an audit trail.
        (output_dir / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    scores = output_dir / "scores.csv"
    write_rows(scores, all_scores)
    provenance["scores_sha256"] = sha256(scores)
    (output_dir / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    return scores


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bagel-root", type=Path)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--scores", type=Path, help="Regenerate report from recorded scores without rerunning fits")
    parser.add_argument("--seed", type=int, default=1729)
    parser.add_argument("--permutation-round", type=int, default=2)
    parser.add_argument("--bootstraps", type=int, default=500)
    args = parser.parse_args()
    if not 0 <= args.seed < 2**32 or min(args.permutation_round, args.bootstraps) < 1:
        parser.error("seed must fit uint32 and rounds/bootstraps must be positive")
    if args.scores is None and args.bagel_root is None:
        parser.error("bagel-root is required unless scores are supplied")
    args.out_dir = args.out_dir.resolve()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["OMP_NUM_THREADS"] = "1"
    scores = args.scores or benchmark(args.bagel_root.resolve(), args.out_dir,
                                     seed=args.seed, rounds=args.permutation_round)
    for row in report(scores, args.out_dir, seed=args.seed, bootstraps=args.bootstraps):
        print(f"{row['dataset']} {row['method']}: AP={row['average_precision']:.4f}", flush=True)


if __name__ == "__main__":
    main()
