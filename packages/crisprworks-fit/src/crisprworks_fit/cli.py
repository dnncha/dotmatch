"""Run the upstream MLE workflow with an optional, measured accelerator."""

import argparse
from contextlib import nullcontext
import hashlib
import json
from pathlib import Path
import platform
import sys
import time

from . import __version__


def positive_int(value):
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be at least 1")
    return parsed


def seed_value(value):
    parsed = int(value)
    if not 0 <= parsed < 2**32:
        raise argparse.ArgumentTypeError("must be between 0 and 4294967295")
    return parsed


def parser():
    result = argparse.ArgumentParser(
        prog="crisprworks-fit", description="CRISPRWorks Fit: experimental MAGeCK2 MLE acceleration",
    )
    result.add_argument("--version", action="version", version=f"CRISPRWorks Fit {__version__}")
    subcommands = result.add_subparsers(dest="subcmd", required=True)
    from mageck2.argsParser import arg_mle
    arg_mle(subcommands)
    mle = subcommands.choices["mle"]
    mle.add_argument("--backend", choices=("accelerated", "reference"), default="accelerated")
    mle.add_argument("--kernel", choices=("auto", "native", "numpy"), default="auto")
    mle.add_argument("--seed", type=seed_value, default=0, help="Permutation RNG seed; default 0")
    mle.add_argument("--blas-threads", type=positive_int, default=1,
                     help="BLAS threads per process; default 1")
    mle.add_argument("--write-fit-details", action="store_true",
                     help="Write full-precision estimates and p-values for numerical comparison")
    demo = subcommands.add_parser("demo", help="Run a reproducible synthetic screen")
    demo.add_argument("--out-dir", type=Path, required=True)
    demo.add_argument("--backend", choices=("accelerated", "reference"), default="accelerated")
    demo.add_argument("--kernel", choices=("auto", "native", "numpy"), default="auto")
    demo.add_argument("--threads", type=positive_int, default=1)
    return result


def synthetic_screen(directory, *, genes=24, guides=6, seed=1729):
    """Generate NB counts for smoke checks; not a biological validation dataset."""
    import numpy as np
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    counts_path = directory / "counts.tsv"
    design_path = directory / "design.tsv"
    generator = np.random.default_rng(seed)
    samples = ("baseline1", "baseline2", "control1", "control2", "treated1", "treated2")
    design_path.write_text("Samples\tbaseline\ttreatment\n" + "".join(
        f"{sample}\t1\t{int(i >= 4)}\n" for i, sample in enumerate(samples)
    ))
    with counts_path.open("w") as stream:
        stream.write("sgRNA\tGene\t" + "\t".join(samples) + "\n")
        for gene in range(genes):
            effect = -1.1 if gene % 6 == 0 else (0.6 if gene % 6 == 1 else 0.0)
            for guide in range(guides):
                mean = generator.uniform(100, 1000)
                values = []
                for sample in range(len(samples)):
                    expected = mean * np.exp(effect if sample >= 4 else 0)
                    values.append(str(generator.negative_binomial(25, 25 / (25 + expected))))
                stream.write(f"g{gene}_{guide}\tGENE{gene}\t" + "\t".join(values) + "\n")
    return counts_path, design_path


def input_record(value):
    if value is None:
        return None
    path = Path(value)
    try:
        is_file = path.is_file()
    except OSError:
        is_file = False
    if is_file:
        digest = hashlib.sha256()
        with path.open("rb") as stream:
            for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(chunk)
        return {"path": str(path.resolve()), "sha256": digest.hexdigest()}
    return {"inline": value, "sha256": hashlib.sha256(value.encode()).hexdigest()}


def run(args):
    import numpy as np
    import scipy
    import mageck2
    from mageck2.mlemageck import mageckmle_main
    from mageck2.mledesignmat import DesignMatCache
    from threadpoolctl import threadpool_limits, threadpool_info
    from .backend import accelerated, verify_upstream, REFERENCE_COMMIT

    verify_upstream()
    from .kernels import resolved_engine
    kernel = resolved_engine(args.kernel) if args.backend == "accelerated" else "reference"
    if args.threads < 1 or args.permutation_round < 1:
        raise ValueError("--threads and --permutation-round must be at least 1")
    if args.max_sgrnapergene_permutation < 2:
        raise ValueError("--max-sgrnapergene-permutation must be at least 2")
    # These workflow branches need separate scientific fixtures before we
    # advertise them. Refuse the options rather than implying validation.
    if args.cnv_norm is not None or args.cnv_est is not None or args.debug_gene is not None:
        raise ValueError("CNV correction and --debug-gene are not supported in this prototype")
    prefix = Path(args.output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    manifest_path = Path(str(prefix) + ".crisprworks.json")
    manifest = {
        "schema_version": 1, "status": "running", "tool": "CRISPRWorks Fit",
        "version": __version__, "backend": args.backend, "kernel": kernel,
        "mageck2_version": mageck2.__version__, "reference_commit": REFERENCE_COMMIT,
        "python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__,
        "platform": platform.platform(), "options": vars(args).copy(),
        "inputs": {name: input_record(getattr(args, name)) for name in (
            "count_table", "design_matrix", "sgrna_efficiency", "control_sgrna", "control_gene",
        )},
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    rng_state = np.random.get_state()
    cache = DesignMatCache.cache
    started = time.perf_counter()
    try:
        np.random.seed(args.seed)
        DesignMatCache.cache = {}
        context = accelerated(args.blas_threads, args.kernel) if args.backend == "accelerated" else nullcontext()
        with threadpool_limits(limits=args.blas_threads), context:
            manifest["blas"] = threadpool_info()
            result = mageckmle_main(parsedargs=args)
        manifest["status"] = "complete"
        manifest["outputs"] = {
            name: input_record(str(prefix) + suffix) for name, suffix in (
                ("gene_summary", ".gene_summary.txt"), ("sgrna_summary", ".sgrna_summary.txt"),
            )
        }
        if args.write_fit_details:
            details_path = Path(str(prefix) + ".fit-details.json")
            fields = (
                "beta_estimate", "beta_zscore", "w_estimate", "beta_pval", "beta_pval_fdr",
                "beta_permute_pval", "beta_permute_pval_fdr", "beta_permute_pval_neg",
                "beta_permute_pval_pos", "beta_permute_pval_neg_fdr", "beta_permute_pval_pos_fdr",
            )
            details = {"schema_version": 1, "genes": {
                name: {"guides": gene.nb_count.shape[1], **{
                    field: np.asarray(getattr(gene, field)).tolist() for field in fields
                }} for name, gene in result[0].items()
            }}
            details_path.write_text(json.dumps(details, indent=2, allow_nan=False) + "\n")
            manifest["outputs"]["fit_details"] = input_record(str(details_path))
        return result
    except BaseException as error:
        manifest["status"] = "failed"
        manifest["error"] = f"{type(error).__name__}: {error}"
        raise
    finally:
        np.random.set_state(rng_state)
        DesignMatCache.cache = cache
        manifest["elapsed_seconds"] = time.perf_counter() - started
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")


def main(argv=None):
    arguments = sys.argv[1:] if argv is None else argv
    cli = parser()
    args = cli.parse_args(arguments)
    if args.subcmd == "demo":
        counts, design = synthetic_screen(args.out_dir)
        args = cli.parse_args([
            "mle", "-k", str(counts), "-d", str(design), "-n", str(args.out_dir / "screen"),
            "--backend", args.backend, "--kernel", args.kernel, "--threads", str(args.threads),
        ])
    try:
        run(args)
    except (ValueError, RuntimeError) as error:
        cli.exit(1, f"crisprworks-fit: {error}\n")
    print(f"CRISPRWorks Fit wrote {args.output_prefix}.gene_summary.txt")
    return 0
