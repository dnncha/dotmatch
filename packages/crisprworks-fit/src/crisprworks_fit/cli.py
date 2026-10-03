"""Run experimental CRISPR gene-effect workflows and reproducible demos."""

import argparse
from contextlib import nullcontext
import hashlib
import json
import math
from pathlib import Path
import platform
import sys
import time
import tempfile

from . import __version__


def positive_int(value):
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be at least 1")
    return parsed


def positive_float(value):
    parsed = float(value)
    if not math.isfinite(parsed) or parsed <= 0:
        raise argparse.ArgumentTypeError("must be finite and positive")
    return parsed


def seed_value(value):
    parsed = int(value)
    if not 0 <= parsed < 2**32:
        raise argparse.ArgumentTypeError("must be between 0 and 4294967295")
    return parsed


def parser():
    result = argparse.ArgumentParser(
        prog="crisprworks-fit", description="CRISPRWorks Fit: experimental gene-level CRISPR inference",
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
    mle.add_argument("--permutation-pvalues", choices=("legacy", "finite"), default="legacy",
                     help="legacy preserves MAGeCK tails; finite counts ties and uses add-one Monte Carlo tails")
    demo = subcommands.add_parser("demo", help="Run a reproducible synthetic screen")
    demo.add_argument("--out-dir", type=Path, required=True)
    demo.add_argument("--backend", choices=("accelerated", "reference"), default="accelerated")
    demo.add_argument("--kernel", choices=("auto", "native", "numpy"), default="auto")
    demo.add_argument("--threads", type=positive_int, default=1)
    demo.add_argument("--model", choices=("mle", "joint"), default="mle",
                      help="mle fits the original demo; joint also calibrates whole control genes")
    calibration = subcommands.add_parser("calibrate", help="Calibrate frozen effects against whole negative-control genes")
    calibration.add_argument("--fit-details", type=Path, required=True)
    calibration.add_argument("--control-gene", type=Path, required=True)
    calibration.add_argument("-n", "--output-prefix", type=Path, required=True)
    calibration.add_argument("--score", choices=("effect", "z"), default="effect")
    calibration.add_argument("--adjust", choices=("bh", "by"), default="by")
    calibration.add_argument("--family", choices=("global", "condition"), default="global")
    calibration.add_argument("--fdr-alpha", type=positive_float, default=.05,
                             help="Resolution diagnostic threshold in (0, 1]; does not change p/q values")
    calibration.add_argument("--min-controls", type=positive_int, default=20)
    calibration.add_argument("--max-guides", type=positive_int, default=40,
                             help="Exclusive MLE permutation guide-count limit; does not exclude joint fits")
    joint = subcommands.add_parser("joint", help="Learn guide efficacy jointly across replicated screens (JACKS-derived)")
    joint.add_argument("-k", "--counts", type=Path, required=True)
    joint.add_argument("--sample-map", type=Path, required=True)
    joint.add_argument("--control-gene", type=Path)
    joint.add_argument("-n", "--output-prefix", type=Path, required=True)
    joint.add_argument("--max-iterations", type=positive_int, default=50)
    joint.add_argument("--tolerance", type=positive_float, default=.1,
                       help="Absolute change in the upstream stopping statistic; default 0.1")
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


def validate_mle_counts(path):
    """Validate the same literal TSV/CSV fields consumed by pinned MAGeCK2."""
    delimiter = ',' if str(path).upper().endswith('.CSV') else '\t'
    seen = set()
    with Path(path).open(encoding='utf-8') as stream:
        header = next(stream, '').strip().split(delimiter)
        if len(header) < 3 or len(set(header)) != len(header) or any(not value.strip() for value in header):
            raise ValueError('MLE count table requires distinct nonempty guide, gene and sample columns')
        for line_number, line in enumerate(stream, 2):
            row = line.strip().split(delimiter)
            if len(row) != len(header) or not row[0].strip() or not row[1].strip():
                raise ValueError(f'Malformed MLE count row {line_number}')
            if row[0] in seen:
                raise ValueError(f'Duplicate guide ID on MLE count row {line_number}')
            seen.add(row[0])
            try:
                valid = all(math.isfinite(float(value)) and float(value) >= 0 for value in row[2:])
            except ValueError:
                valid = False
            if not valid:
                raise ValueError(f'MLE counts must be finite and nonnegative (row {line_number})')
    if not seen:
        raise ValueError('MLE count table requires at least one guide')


def run(args):
    import numpy as np
    import scipy
    import mageck2
    from mageck2.mlemageck import mageckmle_main
    from mageck2.mledesignmat import DesignMatCache
    from threadpoolctl import threadpool_limits, threadpool_info
    from .backend import accelerated, permutation_calibration, verify_upstream, REFERENCE_COMMIT

    verify_upstream()
    from .kernels import resolved_engine
    kernel = resolved_engine(args.kernel) if args.backend == "accelerated" else "reference"
    if args.threads < 1 or args.permutation_round < 1:
        raise ValueError("--threads and --permutation-round must be at least 1")
    if args.max_sgrnapergene_permutation < 2:
        raise ValueError("--max-sgrnapergene-permutation must be at least 2")
    if (args.permutation_pvalues == "finite" and args.no_permutation_by_group
            and (args.control_gene is not None or args.control_sgrna is not None)):
        raise ValueError("Finite control-guide inference requires grouped permutations; remove --no-permutation-by-group")
    # These workflow branches need separate scientific fixtures before we
    # advertise them. Refuse the options rather than implying validation.
    if args.cnv_norm is not None or args.cnv_est is not None or args.debug_gene is not None:
        raise ValueError("CNV correction and --debug-gene are not supported in this prototype")
    prefix = Path(args.output_prefix)
    # Resolve aliases before any output is created, including the running manifest.
    inputs = set()
    for name in ("count_table", "design_matrix", "sgrna_efficiency", "control_sgrna", "control_gene"):
        value = getattr(args, name, None)
        if isinstance(value, (str, Path)):
            candidate = Path(value)
            try:
                if candidate.is_file():
                    inputs.add(candidate.resolve())
            except OSError:
                pass
    suffixes = (".crisprworks.json", ".gene_summary.txt", ".sgrna_summary.txt", ".fit-details.json")
    from .publication import reject_input_collisions, publish_file_bundle
    reject_input_collisions((Path(str(prefix) + suffix) for suffix in suffixes), inputs)
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
    calibration_groups = []
    if args.permutation_pvalues == "finite":
        manifest["permutation_calibration"] = {
            "method": "inclusive-add-one", "two_sided": "min(1, 2 * min(upper, lower))",
            "invalid_fit_policy": "all tails equal one", "skipped_fit_policy": "all tails equal one",
            "null_source": "control guides" if (args.control_sgrna is not None or args.control_gene is not None) else "all guides",
            "scope": "finite-tail safeguard; pooled-guide exchangeability and biological FDR are not established",
            "groups": calibration_groups,
        }
    rng_state = np.random.get_state()
    cache = DesignMatCache.cache
    started = time.perf_counter()
    staging = tempfile.TemporaryDirectory(prefix=f".{prefix.name}.mle-", dir=prefix.parent)
    staged_prefix = Path(staging.name) / "result"
    original_prefix = args.output_prefix
    original_inputs = {name: getattr(args, name) for name in manifest["inputs"]}
    try:
        # Upstream reads immutable copies; provenance hashes exactly those bytes.
        for name, record in manifest["inputs"].items():
            if isinstance(record, dict) and "path" in record:
                payload = Path(original_inputs[name]).read_bytes()
                # MAGeCK selects CSV delimiters from the filename suffix.
                snapshot = Path(staging.name) / (name + Path(original_inputs[name]).suffix)
                snapshot.write_bytes(payload)
                record["sha256"] = hashlib.sha256(payload).hexdigest()
                setattr(args, name, str(snapshot))
        args.output_prefix = str(staged_prefix)
        validate_mle_counts(args.count_table)
        np.random.seed(args.seed)
        DesignMatCache.cache = {}
        context = accelerated(args.blas_threads, args.kernel) if args.backend == "accelerated" else nullcontext()
        with threadpool_limits(limits=args.blas_threads), context, permutation_calibration(
            args.permutation_pvalues, max_guides=args.max_sgrnapergene_permutation,
            diagnostics=calibration_groups,
        ):
            manifest["blas"] = threadpool_info()
            result = mageckmle_main(parsedargs=args)
        manifest["status"] = "complete"
        outputs = {name: (Path(str(staged_prefix) + suffix), Path(str(prefix) + suffix))
                   for name, suffix in (("gene_summary", ".gene_summary.txt"),
                                        ("sgrna_summary", ".sgrna_summary.txt"))}
        if args.write_fit_details:
            details_path = Path(str(staged_prefix) + ".fit-details.json")
            fields = (
                "beta_estimate", "beta_zscore", "w_estimate", "beta_pval", "beta_pval_fdr",
                "beta_permute_pval", "beta_permute_pval_fdr", "beta_permute_pval_neg",
                "beta_permute_pval_pos", "beta_permute_pval_neg_fdr", "beta_permute_pval_pos_fdr",
            )
            n_conditions = args.design_matrix.shape[1] - 1
            condition_labels = (list(args.beta_labels[1:]) if args.beta_labels is not None else
                                [f"beta_{i + 1}" for i in range(n_conditions)])
            details = {"schema_version": 1, "model": "MAGeCK2 MLE",
                       "effect_units": "natural-log beta coefficient", "conditions": condition_labels,
                       "genes": {
                name: {"guides": gene.nb_count.shape[1], **{
                    field: np.asarray(getattr(gene, field)).tolist() for field in fields
                }} for name, gene in result[0].items()
            }}
            details_path.write_text(json.dumps(details, indent=2, allow_nan=False) + "\n")
            outputs["fit_details"] = (details_path, Path(str(prefix) + ".fit-details.json"))
        manifest["elapsed_seconds"] = time.perf_counter() - started
        remove = () if args.write_fit_details else (Path(str(prefix) + ".fit-details.json"),)
        publish_file_bundle(outputs, manifest_path, manifest, remove)
        return result
    except BaseException as error:
        manifest["status"] = "failed"
        manifest["error"] = f"{type(error).__name__}: {error}"
        raise
    finally:
        np.random.set_state(rng_state)
        DesignMatCache.cache = cache
        manifest["elapsed_seconds"] = time.perf_counter() - started
        args.output_prefix = original_prefix
        for name, value in original_inputs.items():
            setattr(args, name, value)
        staging.cleanup()


def main(argv=None):
    arguments = sys.argv[1:] if argv is None else argv
    cli = parser()
    args = cli.parse_args(arguments)
    if args.subcmd in ("calibrate", "joint"):
        if args.subcmd == "calibrate":
            from .calibration import run_calibration as operation
        else:
            from .joint import run_joint as operation
        try:
            output = operation(args)
        except (ValueError, KeyError, OSError, TypeError) as error:
            cli.exit(1, f"crisprworks-fit: {error}\n")
        print(f"CRISPRWorks Fit wrote {output}")
        return 0
    if args.subcmd == "demo":
        if args.model == "joint":
            if args.backend != "accelerated" or args.kernel != "auto" or args.threads != 1:
                cli.error("--backend, --kernel and --threads select MLE demo behavior; omit them with --model joint")
            from .demo import run_joint_demo
            try:
                output = run_joint_demo(args.out_dir)
            except (ValueError, OSError, KeyError, TypeError) as error:
                cli.exit(1, f"crisprworks-fit: {error}\n")
            print(f"CRISPRWorks Fit wrote {output}; see {args.out_dir / 'screen.calibration.tsv'} for control-calibrated results")
            return 0
        counts, design = synthetic_screen(args.out_dir)
        args = cli.parse_args([
            "mle", "-k", str(counts), "-d", str(design), "-n", str(args.out_dir / "screen"),
            "--backend", args.backend, "--kernel", args.kernel, "--threads", str(args.threads),
        ])
    try:
        run(args)
    except (ValueError, RuntimeError, OSError) as error:
        cli.exit(1, f"crisprworks-fit: {error}\n")
    print(f"CRISPRWorks Fit wrote {args.output_prefix}.gene_summary.txt")
    return 0
