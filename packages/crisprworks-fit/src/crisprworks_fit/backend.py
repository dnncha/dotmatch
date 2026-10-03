"""Scoped acceleration of the tested MAGeCK2 implementation."""

from concurrent.futures import ProcessPoolExecutor
from contextlib import contextmanager
import hashlib
import inspect
import multiprocessing
import threading

import numpy as np

from . import kernels

REFERENCE_COMMIT = "630aea0b6fc152a81006435f21d629297275c911"
FUNCTION_HASHES = {
    "em_whileloop": "f744a81f54461a3d6b06799f50a9bb74398ddce265d6141fcf5f6c2bb31e5896",
    "iteratenbem": "8c04d689c811a8829e6586b6f25b86ee9e6ce3bf944715517a22866b646f075e",
    "getloglikelihood2": "b75f0f86f4926737b3fde9646082c663dd9e26b01b9fa2a3d2fa9cca83bed88d",
    "runem_multiproc": "c07a3b55b2dec0b3293b7bdc6409ec526e14487b0698a7b3ee10f6d2fc870703",
    "assign_p_value_from_permuted_beta": "3aa3deb77c5d22389a0d09adabcd0e5214422790e63e5e6619b7b20a4484715c",
}
_lock = threading.RLock()
_active = False
_calibration_active = False


def verify_upstream():
    import mageck2
    from mageck2 import mleem, mlemultiprocessing

    if mageck2.__version__ != "0.3.0":
        raise RuntimeError("CRISPRWorks Fit currently requires mageck2==0.3.0")
    for module, names in (
        (mleem, ("em_whileloop", "iteratenbem", "getloglikelihood2")),
        (mlemultiprocessing, ("runem_multiproc", "assign_p_value_from_permuted_beta")),
    ):
        for name in names:
            digest = hashlib.sha256(inspect.getsource(getattr(module, name)).encode()).hexdigest()
            if digest != FUNCTION_HASHES[name]:
                raise RuntimeError(f"Untested MAGeCK2 implementation of {name}; use the pinned release")


@contextmanager
def permutation_calibration(mode="legacy", *, max_guides=None, diagnostics=None):
    """Select inference after entering the reference or accelerated context.

    Like acceleration, this patch is scoped and process-global. The reference
    mode is still MAGeCK fitting when finite tails are explicitly requested;
    its inference then deliberately differs from unmodified MAGeCK.
    """
    global _calibration_active
    if mode not in ("legacy", "finite"):
        raise ValueError("Permutation p-values must be legacy or finite")
    if mode == "legacy":
        yield
        return
    with _lock:
        if _calibration_active:
            raise RuntimeError("Nested permutation calibration contexts are unsupported")
        if not _active:
            verify_upstream()
        from mageck2 import mlemultiprocessing
        from .inference import assign_finite_permutation_pvalues
        original = mlemultiprocessing.assign_p_value_from_permuted_beta
        def assign(null, genes):
            return assign_finite_permutation_pvalues(
                null, genes, max_guides=max_guides, diagnostics=diagnostics,
            )
        _calibration_active = True
        mlemultiprocessing.assign_p_value_from_permuted_beta = assign
        try:
            yield
        finally:
            mlemultiprocessing.assign_p_value_from_permuted_beta = original
            _calibration_active = False


def _worker_init(blas_threads, kernel):
    verify_upstream()
    kernels._engine = kernels.resolved_engine(kernel)
    from mageck2 import mleem
    from mageck2.mledesignmat import DesignMatCache
    DesignMatCache.cache = {}
    mleem.em_whileloop = kernels.em_whileloop
    from threadpoolctl import threadpool_limits
    threadpool_limits(limits=blas_threads)


def _fit_one(task):
    from mageck2.mleem import iteratenbem
    from mageck2.mledesignmat import DesignMatCache
    gene_id, gene, options = task
    # Upstream's cache key only contains guide count; explicitly refresh the
    # design so repeated library calls with different designs remain safe.
    DesignMatCache.save_record(gene.design_mat, gene.nb_count.shape[1])
    iteratenbem(gene, **options)
    return gene_id, gene


class Runner:
    """Reuse a worker pool across fitting stages and permutation rounds."""

    def __init__(self, blas_threads=1, kernel="auto"):
        self.executor = None
        self.workers = None
        self.blas_threads = blas_threads
        self.kernel = kernels.resolved_engine(kernel)

    def close(self):
        if self.executor is not None:
            self.executor.shutdown(wait=True, cancel_futures=True)
            self.executor = None

    def __call__(self, allgenedict, args, nproc=1, argsdict=None):
        if nproc < 1:
            raise ValueError("--threads must be at least 1")
        options = {} if argsdict is None else argsdict
        tasks = []
        for gene_id, gene in allgenedict.items():
            n_guides = gene.nb_count.shape[1]
            n_conditions = gene.design_mat.shape[1] - 1
            # Preserve upstream's inclusive skip threshold.
            if n_guides >= args.max_sgrnapergene_permutation:
                gene.beta_estimate = np.zeros(n_guides + n_conditions)
                gene.beta_zscore = np.zeros(n_conditions)
                for name in ("beta_pval", "beta_pval_pos", "beta_pval_neg"):
                    setattr(gene, name, np.ones(n_conditions))
            else:
                tasks.append((gene_id, gene, options))
        if nproc == 1:
            for task in tasks:
                gene_id, fitted = _fit_one(task)
                allgenedict[gene_id] = fitted
            return
        if self.executor is None or self.workers != nproc:
            self.close()
            self.workers = nproc
            self.executor = ProcessPoolExecutor(
                max_workers=nproc, mp_context=multiprocessing.get_context("spawn"),
                initializer=_worker_init,
                initargs=(self.blas_threads, self.kernel),
            )
        # map yields in input order, preserving gene order for permutations.
        chunk = max(1, len(tasks) // (nproc * 4))
        for gene_id, fitted in self.executor.map(_fit_one, tasks, chunksize=chunk):
            allgenedict[gene_id] = fitted


@contextmanager
def accelerated(blas_threads=1, kernel="auto"):
    """Enable kernels for one workflow, restoring upstream state on exit.

    Patches are process-global: concurrent calls to unwrapped MAGeCK2 in the
    same interpreter are unsupported. Use separate processes for other work.
    """
    global _active
    with _lock:
        if _active:
            raise RuntimeError("Nested accelerator contexts are unsupported")
        verify_upstream()
        from mageck2 import mleem, mlemultiprocessing
        from mageck2.mledesignmat import DesignMatCache
        original = (mleem.em_whileloop, mlemultiprocessing.runem_multiproc,
                    mlemultiprocessing.assign_p_value_from_permuted_beta, DesignMatCache.cache)
        runner = Runner(blas_threads, kernel)
        previous_engine = kernels._engine
        kernels._engine = runner.kernel
        _active = True
        DesignMatCache.cache = {}
        mleem.em_whileloop = kernels.em_whileloop
        mlemultiprocessing.runem_multiproc = runner
        mlemultiprocessing.assign_p_value_from_permuted_beta = kernels.assign_p_value_from_permuted_beta
        try:
            yield
        finally:
            try:
                runner.close()
            finally:
                (mleem.em_whileloop, mlemultiprocessing.runem_multiproc,
                 mlemultiprocessing.assign_p_value_from_permuted_beta, DesignMatCache.cache) = original
                kernels._engine = previous_engine
                _active = False
