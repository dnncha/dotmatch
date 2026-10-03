"""Plot annotated-class call rates; these are not genome-wide FDR estimates."""

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metrics", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True, help="SVG output path; a PDF is also written")
    args = parser.parse_args()
    with args.metrics.open() as stream:
        rows = {(row["dataset"], row["method"]): row for row in csv.DictReader(stream)}
    datasets = ("hap1-tko", "hap1-tkov3")
    lanes = (("fit_finite_depletion_p", "All-guide null", "#737d88"),
             ("fit_control_depletion_p", "Held-out control null", "#246b9e"))
    figure, axes = plt.subplots(1, 2, figsize=(9, 3.7), layout="constrained")
    for axis, field, title in zip(axes, ("essential_call_rate_at_nominal_fdr_005", "nonessential_call_rate_at_nominal_fdr_005"),
                                 ("Reference essential sensitivity", "Reference nonessential call rate")):
        for i, (method, label, color) in enumerate(lanes):
            values = [float(rows[dataset, method][field]) for dataset in datasets]
            positions = np.arange(2) + (i - .5) * .34
            axis.bar(positions, values, width=.3, color=color, label=label)
            for position, value in zip(positions, values):
                axis.annotate(f"{100 * value:.1f}%", (position, value), xytext=(0, 4),
                              textcoords="offset points", ha="center", fontsize=9)
        axis.set_xticks(range(2), ("HAP1 TKO", "HAP1 TKOv3"))
        axis.set_ylim(0, 1.12 if "essential_call" in field and "nonessential" not in field else .12)
        axis.set_title(title, fontsize=11)
        axis.set_ylabel("Fraction of annotated class called")
        axis.spines[["top", "right"]].set_visible(False)
    handles, labels = axes[0].get_legend_handles_labels()
    figure.legend(handles, labels, frameon=False, loc="lower center", ncol=2,
                  bbox_to_anchor=(.5, -.07), fontsize=10)
    figure.suptitle("Fit: negative effects with directional FDR ≤ 0.05\nFive external gene folds; complete four-guide cohorts", fontsize=12)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.out, bbox_inches="tight")
    figure.savefig(args.out.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(figure)


if __name__ == "__main__":
    main()
