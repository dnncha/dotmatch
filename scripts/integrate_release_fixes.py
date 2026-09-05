#!/usr/bin/env python3
"""One-shot, hash-guarded integration; removed from the final release tree."""
import hashlib
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
EXPECTED = {
    'python/dotmatch/core.py': 'edad2877278ae0ceb9bcfe0143b3bd81a6b96709',
    'python/dotmatch/cli.py': 'c66c3f3dfb66dfadb5dd79ab35c37efc23a1d00c',
    'python/dotmatch/tl.py': '9cc688db3b1f5ea10ba02b04711d14e28d439cbe',
    'README.md': 'e4f1c3fbc82d30935bbc74ebc639af134eaec7cc',
    'scripts/check_site_browser.py': '8f1e014e674163b4dff7a13add500bce80781dff',
}
for name, sha in EXPECTED.items():
    data = (ROOT / name).read_bytes()
    assert hashlib.sha1(f'blob {len(data)}\0'.encode() + data).hexdigest() == sha, f'Concurrent edit: {name}'

def replace(text, old, new, count=1):
    assert text.count(old) == count, f'Expected {count} occurrence(s): {old[:100]!r}'
    return text.replace(old, new)

def edit(path, function):
    file = ROOT / path
    file.write_text(function(file.read_text(encoding='utf-8')), encoding='utf-8')

def core(text):
    for name, alias, flag in [('pandas', 'pd', 'PANDAS'), ('polars', 'pl', 'POLARS'), ('anndata', 'ad', 'ANNDATA')]:
        pattern = rf'^try:\n    import {name} as {alias}[^\n]*\n    _HAS_{flag} = True\nexcept Exception:[^\n]*\n    {alias} = None[^\n]*\n    _HAS_{flag} = False'
        text, count = re.subn(pattern, f'{alias}, _HAS_{flag} = optional_module("{name}")', text, flags=re.M)
        assert count == 1, name
    text = replace(text, 'from typing import Any, Iterable, Iterator, Sequence, TextIO\n', 'from typing import Any, Iterable, Iterator, Sequence, TextIO\n\nfrom ._optional import optional_module\n')
    text = replace(text, 'return gzip.open(path, mode)', 'return gzip.open(path, mode, encoding="utf-8", newline="")')
    text = replace(text, 'def iter_fastq(path: str | Path) -> Iterator[FastqRecord]:\n    """Yield FASTQ records from plain or gzipped FASTQ."""', '''def iter_fastq(path: str | Path, *, content_digest: Any | None = None) -> Iterator[FastqRecord]:
    """Yield FASTQ records, optionally hashing original decompressed UTF-8 bytes.

    The optional hashlib-compatible object is updated before normalization.
    Consume the complete iterator before treating its digest as a full-input
    checksum; this is not the checksum of a compressed file.
    """''')
    text = replace(text, '            qual = fh.readline()\n            if not seq or not plus or not qual:', '''            qual = fh.readline()
            if content_digest is not None:
                for line in (header, seq, plus, qual):
                    content_digest.update(line.encode("utf-8"))
            if not seq or not plus or not qual:''')
    text = replace(text, '            yield FastqRecord(header[1:].split()[0], seq, qual)', '''            identifiers = header[1:].split()
            if not identifiers:
                raise ValueError(f"invalid FASTQ record in {path}: missing read identifier")
            yield FastqRecord(identifiers[0], seq, qual)''')
    text = replace(text, '''    if _HAS_PANDAS and hasattr(targets, "columns"):
        return targets_from_dataframe(targets)
    if _HAS_POLARS and pl is not None and isinstance(targets, pl.DataFrame):
        return targets_from_dataframe(targets)''', '''    # Do not import an optional dataframe stack to inspect ordinary lists.
    if hasattr(targets, "columns"):
        return targets_from_dataframe(targets)''')
    text = replace(text, '    seqs = data[seq_col].astype(str).tolist()', '    seqs = data[seq_col].astype(str).str.upper().tolist()')
    text = replace(text, '''        if target_names is not None and 0 <= r.target_index < len(target_names):
            row["target_name"] = target_names[r.target_index]''', '''        if target_names is not None:
            # A candidate index is diagnostic information, not an assignment.
            row["target_name"] = (
                target_names[r.target_index]
                if r.status == MATCH_UNIQUE and 0 <= r.target_index < len(target_names)
                else ""
            )''')
    return text

edit('python/dotmatch/core.py', core)
edit('python/dotmatch/tl.py', lambda text: replace(replace(text, '''try:
    import anndata as ad
    import pandas as pd
    _HAS_ANNDATA = True
except Exception:  # noqa: BLE001
    ad = None
    pd = None
    _HAS_ANNDATA = False''', '''from ._optional import optional_module

ad, _HAS_ANNDATA = optional_module("anndata")
pd, _HAS_PANDAS = optional_module("pandas")'''), '    if not _HAS_ANNDATA:', '    if not (_HAS_ANNDATA and _HAS_PANDAS):'))
edit('python/dotmatch/cli.py', lambda text: replace(text, '''                    "Draft project created; review inference_report.json, then rerun with "
                    "--accept-inference or use dotmatch assay start after setting status = \\"ready\\"."'''.replace('\\\\', '\\'), '''                    "Draft project created; review inference_report.json and assay.toml. "
                    "After confirming the settings, change status = \\"draft\\" to \\"ready\\" "
                    "in assay.toml, then run dotmatch assay start with that file."'''.replace('\\\\', '\\')))

def readme(text):
    text = replace(text, '''Keep downstream screen statistics in the workflow you already use. DotMatch is
an assignment and counting tool, not a genome aligner, basecaller, cell/UMI
pipeline or gene-level hit-calling package.''', '''Keep downstream screen statistics in the workflow you already use. DotMatch is not a genome aligner,
basecaller, cell/UMI pipeline or gene-level hit-calling package.''')
    text = replace(text, '''Bioconda and its generated BioContainers images can lag PyPI/GHCR. Check the
installed version.''', '''Bioconda and its generated BioContainers images can lag PyPI/GHCR. When a
newly tagged version has not reached Bioconda yet, use PyPI or the source build.
Check the installed version. Review the [packaging details](https://dotmatch.readthedocs.io/en/latest/packaging.html)
for platform and container verification.''')
    text = replace(text, '''Review `crispr_screen/inference_report.json` and `assay.toml`: confirm the guide
window, orientation, library and sample files. Then run and review:''', '''This creates a draft project. Review `crispr_screen/inference_report.json` and
`assay.toml`: confirm the guide window, orientation, library and sample files.
After confirming the settings, change the top-level `status = "draft"` to
`status = "ready"` in `assay.toml`, then run and review:''')
    text = replace(text, '## Reproduce the evidence\n', '''## Reproduce the evidence

The [benchmark reports](https://dotmatch.readthedocs.io/en/latest/benchmarks/README.html)
include commands, hardware and assignment rules. Those reports cover the tested workloads;
they are not universal speed or biological-accuracy guarantees.
''')
    return text
edit('README.md', readme)
edit('scripts/check_site_browser.py', lambda text: replace(replace(text, 'from playwright.sync_api import sync_playwright', 'from playwright.sync_api import expect, sync_playwright'), 'assert "A, C, G and T" in page.get_by_role("alert").inner_text()', 'expect(page.locator("#library-error")).to_contain_text("A, C, G and T")'))
print('Integrated exact-scope package fixes and CI corrections without replacing newer sensitivity or website work.')
