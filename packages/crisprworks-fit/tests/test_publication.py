"""Input preservation regressions for result publication."""

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from crisprworks_fit.cli import parser, run
from crisprworks_fit.calibration import run_calibration
from crisprworks_fit.publication import unique_staging_path
from crisprworks_fit.publication import reject_input_collisions, publish_table_bundle
from crisprworks_fit.publication import verify_output_record
from crisprworks_fit.publication import read_text_input
from crisprworks_fit.publication import publish_file_bundle
import subprocess
import sys
from queue import Queue
from threading import Thread


@pytest.mark.parametrize("suffix", [".crisprworks.json", ".gene_summary.txt",
                                    ".sgrna_summary.txt", ".fit-details.json"])
def test_mle_rejects_input_output_collision_before_writing(tmp_path, suffix):
    prefix = tmp_path / "run"
    counts = Path(str(prefix) + suffix)
    original = b"sgRNA\tGene\tbaseline\ttreated\na\tA\t10\t5\n"
    counts.write_bytes(original)
    args = parser().parse_args(["mle", "-k", str(counts), "-d", "1,0;1,1",
                               "-n", str(prefix), "--write-fit-details", "--kernel", "numpy"])
    with pytest.raises(ValueError, match="collides"):
        run(args)
    assert counts.read_bytes() == original
    assert list(tmp_path.iterdir()) == [counts]


def test_calibration_preserves_input_named_like_old_staging_file(tmp_path, capsys):
    prefix = tmp_path / "run"
    details = tmp_path / "run.calibration.tsv.tmp"
    original = json.dumps({"schema_version": 1, "genes": {
        name: {"guides": 4, "effect_estimate": [effect], "effect_zscore": [effect]}
        for name, effect in (("a", 0), ("b", 1), ("test", -10))
    }}).encode()
    details.write_bytes(original)
    controls = tmp_path / "controls.txt"
    controls.write_text("a\nb\n")
    output = run_calibration(SimpleNamespace(fit_details=details, control_gene=controls,
        output_prefix=prefix, score="effect", adjust="by", family="global",
        min_controls=2, max_guides=40))
    assert output.is_file()
    assert details.read_bytes() == original
    manifest = json.loads((tmp_path / "run.calibration.json").read_text())
    assert manifest["inputs"]["fit_details"]["path"] == str(details.resolve())
    assert manifest["effect_model"] == "JACKS-derived joint guide-efficacy inference"
    assert manifest["effect_units"] == "log2 relative abundance"
    assert manifest["resolution"]["directional_hypotheses_with_attainable_q"] == 0
    assert "no calls would not establish absence of effects" in capsys.readouterr().err
    assert not list(tmp_path.glob(".run.*.tmp"))


def test_staging_paths_are_exclusive_and_leave_existing_files_untouched(tmp_path):
    destination = tmp_path / "result.tsv"
    old = tmp_path / "result.tsv.tmp"
    old.write_text("input")
    first, second = unique_staging_path(destination), unique_staging_path(destination)
    try:
        assert first != second and first != old and second != old
        assert first.is_file() and second.is_file()
        assert old.read_text() == "input"
    finally:
        first.unlink()
        second.unlink()


def test_mle_rejects_hard_link_to_input(tmp_path):
    counts = tmp_path / "counts.tsv"
    original = b"sgRNA\tGene\tbaseline\ttreated\na\tA\t10\t5\n"
    counts.write_bytes(original)
    output = tmp_path / "run.crisprworks.json"
    output.hardlink_to(counts)
    args = parser().parse_args(["mle", "-k", str(counts), "-d", "1,0;1,1",
                               "-n", str(tmp_path / "run"), "--kernel", "numpy"])
    with pytest.raises(ValueError, match="collides"):
        run(args)
    assert counts.read_bytes() == original
    assert output.read_bytes() == original


def test_collision_check_handles_symlink_and_hard_link(tmp_path):
    source = tmp_path / "input"
    source.write_text("preserve")
    for alias, link in ((tmp_path / "symbolic", "symlink_to"),
                        (tmp_path / "hard", "hardlink_to")):
        getattr(alias, link)(source)
        with pytest.raises(ValueError, match="collides"):
            reject_input_collisions([alias], [source])


def test_manifest_staging_failure_preserves_previous_bundle(tmp_path, monkeypatch):
    table, manifest = tmp_path / "table", tmp_path / "manifest"
    staged = unique_staging_path(table)
    table.write_text("old table")
    manifest.write_text("old manifest")
    staged.write_text("new table")
    write_text = Path.write_text
    def fail_staging(path, *args, **kwargs):
        if path.name.startswith(".manifest."):
            raise OSError("injected staging failure")
        return write_text(path, *args, **kwargs)
    monkeypatch.setattr(Path, "write_text", fail_staging)
    try:
        with pytest.raises(OSError, match="injected"):
            publish_table_bundle(staged, table, manifest, {})
        assert table.read_text() == "old table"
        assert manifest.read_text() == "old manifest"
        assert not list(tmp_path.glob(".manifest.*.tmp"))
    finally:
        staged.unlink(missing_ok=True)


def test_failure_after_table_replacement_leaves_no_completion_manifest(tmp_path, monkeypatch):
    table, manifest = tmp_path / "table", tmp_path / "manifest"
    staged = unique_staging_path(table)
    staged.write_text("new table")
    table.write_text("old table")
    manifest.write_text("old manifest")
    replace = Path.replace
    def fail_manifest(path, target):
        if target == manifest:
            raise OSError("injected publication failure")
        return replace(path, target)
    monkeypatch.setattr(Path, "replace", fail_manifest)
    with pytest.raises(OSError, match="injected"):
        publish_table_bundle(staged, table, manifest, {})
    assert table.read_text() == "new table"
    assert not manifest.exists()
    assert not list(tmp_path.glob(".manifest.*.tmp"))


def test_bundle_records_replaced_output_path_and_matching_hash(tmp_path):
    import hashlib
    table, manifest = tmp_path / "table", tmp_path / "manifest"
    previous_target = tmp_path / "previous-target"
    previous_target.write_text("leave alone")
    table.symlink_to(previous_target)
    staged = unique_staging_path(table)
    staged.write_text("new table")
    publish_table_bundle(staged, table, manifest, {})
    output = json.loads(manifest.read_text())["output"]
    assert output["path"] == str(table)
    assert output["sha256"] == hashlib.sha256(table.read_bytes()).hexdigest()
    assert previous_target.read_text() == "leave alone"


@pytest.mark.parametrize("record", [None, {}, {"path": "", "sha256": "0" * 64},
                                    {"path": "/missing", "sha256": "z" * 64}])
def test_malformed_companion_records_are_rejected(record):
    with pytest.raises(ValueError, match="Malformed companion"):
        verify_output_record({"output": record})


def test_missing_companion_is_rejected_but_legacy_details_are_supported(tmp_path):
    verify_output_record({"genes": {}})
    with pytest.raises(FileNotFoundError):
        verify_output_record({"output": {"path": str(tmp_path / "missing"), "sha256": "0" * 64}})


def test_calibration_protects_referenced_companion_from_output_collision(tmp_path):
    import hashlib
    companion = tmp_path / "run.calibration.tsv"
    original = b"original fit summary\n"
    companion.write_bytes(original)
    details = tmp_path / "details.json"
    details.write_text(json.dumps({"schema_version": 1, "output": {
        "path": str(companion), "sha256": hashlib.sha256(original).hexdigest()},
        "genes": {name: {"guides": 4, "effect_estimate": [effect], "effect_zscore": [effect]}
                  for name, effect in (("a", 0), ("b", 1), ("test", -10))}}))
    controls = tmp_path / "controls"
    controls.write_text("a\nb\n")
    with pytest.raises(ValueError, match="collides"):
        run_calibration(SimpleNamespace(fit_details=details, control_gene=controls,
            output_prefix=tmp_path / "run", score="effect", adjust="by", family="global",
            min_controls=2, max_guides=40))
    assert companion.read_bytes() == original
    assert not (tmp_path / "run.calibration.json").exists()


@pytest.mark.parametrize("iterations,change", [(2, .01), (1, .1), (1, .2)])
def test_calibration_rejects_diagnostics_contradicting_declared_settings(tmp_path, iterations, change):
    details = tmp_path / "details.json"
    details.write_text(json.dumps({"schema_version": 1, "max_iterations": 1,
        "absolute_bound_tolerance": .1,
        "genes": {name: {"guides": 4, "effect_estimate": [effect], "effect_zscore": [effect],
                         "fit_diagnostics": {"converged": True, "iterations": iterations,
                             "termination_reason": "bound_tolerance", "final_bound_change": change}}
                  for name, effect in (("a", 0), ("b", 1), ("test", -10))}}))
    controls = tmp_path / "controls"
    controls.write_text("a\nb\n")
    with pytest.raises(ValueError, match="contradict"):
        run_calibration(SimpleNamespace(fit_details=details, control_gene=controls,
            output_prefix=tmp_path / "run", score="effect", adjust="by", family="global",
            min_controls=2, max_guides=40))
    assert not list(tmp_path.glob("run.*"))


def test_overlapping_publisher_fails_without_replacing_previous_bundle(tmp_path):
    table, manifest = tmp_path / 'table', tmp_path / 'manifest'
    table.write_text('old table')
    manifest.write_text('old manifest')
    staged = unique_staging_path(table)
    staged.write_text('new table')
    code = ("from pathlib import Path; from crisprworks_fit.publication import publication_lock; "
            "import sys\n"
            f"with publication_lock(Path({str(manifest)!r})):\n"
            " print('locked', flush=True)\n"
            " sys.stdin.readline()\n")
    process = subprocess.Popen([sys.executable, '-c', code], stdin=subprocess.PIPE,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    try:
        ready = Queue()
        Thread(target=lambda: ready.put(process.stdout.readline()), daemon=True).start()
        assert ready.get(timeout=10).strip() == 'locked'
        with pytest.raises(ValueError, match='already in progress'):
            publish_table_bundle(staged, table, manifest, {})
        assert table.read_text() == 'old table'
        assert manifest.read_text() == 'old manifest'
        # Abrupt exit must release the OS lock without deleting the lock file.
        process.kill()
        process.wait(timeout=10)
        publish_table_bundle(staged, table, manifest, {})
        verify_output_record(json.loads(manifest.read_text()))
        assert table.read_text() == 'new table'
    finally:
        if process.poll() is None:
            process.kill()
        process.communicate(timeout=10)
        staged.unlink(missing_ok=True)


def test_changed_inputs_are_rejected_before_publication(tmp_path):
    import hashlib
    source = tmp_path / 'input'
    source.write_text('original')
    metadata = {'inputs': {'counts': {'path': str(source),
                'sha256': hashlib.sha256(source.read_bytes()).hexdigest()}}}
    source.write_text('changed')
    table, manifest = tmp_path / 'table', tmp_path / 'manifest'
    table.write_text('old table')
    manifest.write_text('old manifest')
    staged = unique_staging_path(table)
    staged.write_text('new table')
    try:
        with pytest.raises(ValueError, match='Input changed'):
            publish_table_bundle(staged, table, manifest, metadata)
        assert table.read_text() == 'old table'
        assert manifest.read_text() == 'old manifest'
    finally:
        staged.unlink(missing_ok=True)


def test_output_collision_created_during_staging_is_rechecked(tmp_path, monkeypatch):
    import hashlib
    source = tmp_path / 'input'
    source.write_text('preserve')
    table, manifest = tmp_path / 'table', tmp_path / 'manifest'
    staged = unique_staging_path(table)
    staged.write_text('new table')
    metadata = {'inputs': {'counts': {'path': str(source),
                'sha256': hashlib.sha256(source.read_bytes()).hexdigest()}}}
    original_write = Path.write_text
    def create_alias(path, *args, **kwargs):
        result = original_write(path, *args, **kwargs)
        if path.name.startswith('.manifest.'):
            table.hardlink_to(source)
        return result
    monkeypatch.setattr(Path, 'write_text', create_alias)
    try:
        with pytest.raises(ValueError, match='collides'):
            publish_table_bundle(staged, table, manifest, metadata)
        assert source.read_text() == 'preserve'
        assert table.samefile(source)
        assert not manifest.exists()
    finally:
        staged.unlink(missing_ok=True)


def test_publication_lock_cannot_alias_an_input(tmp_path):
    import hashlib
    source = tmp_path / '.manifest.publish.lock'
    source.write_text('preserve input')
    metadata = {'inputs': {'counts': {'path': str(source),
                'sha256': hashlib.sha256(source.read_bytes()).hexdigest()}}}
    table, manifest = tmp_path / 'table', tmp_path / 'manifest'
    staged = unique_staging_path(table)
    staged.write_text('new table')
    try:
        with pytest.raises(ValueError, match='collides'):
            publish_table_bundle(staged, table, manifest, metadata)
        assert source.read_text() == 'preserve input'
        assert not table.exists() and not manifest.exists()
    finally:
        staged.unlink(missing_ok=True)


def test_input_provenance_hashes_exact_parsed_bytes(tmp_path):
    import hashlib
    source = tmp_path / 'input'
    payload = 'α\r\n'.encode('utf-8')
    source.write_bytes(payload)
    text, record = read_text_input(source)
    source.write_text('different')
    assert text == 'α\r\n'
    assert record['sha256'] == hashlib.sha256(payload).hexdigest()


def mle_fixture(tmp_path):
    counts = tmp_path / 'counts.tsv'
    counts.write_text('sgRNA\tGene\tbaseline\ttreated\na\tA\t10\t5\n')
    prefix = tmp_path / 'result'
    previous = {}
    for suffix in ('.gene_summary.txt', '.sgrna_summary.txt', '.crisprworks.json', '.fit-details.json'):
        path = Path(str(prefix) + suffix)
        path.write_text('previous ' + suffix)
        previous[path] = path.read_bytes()
    args = parser().parse_args(['mle', '-k', str(counts), '-d', '1,0;1,1',
                               '-n', str(prefix), '--kernel', 'numpy'])
    return args, counts, previous


def test_mle_fitting_failure_preserves_all_previous_outputs(tmp_path, monkeypatch):
    args, counts, previous = mle_fixture(tmp_path)
    def fail(parsedargs):
        Path(parsedargs.output_prefix + '.gene_summary.txt').write_text('partial new result')
        raise RuntimeError('injected fitting failure')
    monkeypatch.setattr('mageck2.mlemageck.mageckmle_main', fail)
    with pytest.raises(RuntimeError, match='injected'):
        run(args)
    assert all(path.read_bytes() == content for path, content in previous.items())
    assert args.count_table == str(counts)
    assert args.output_prefix == str(tmp_path / 'result')
    assert not list(tmp_path.glob('.result.mle-*'))


def test_successful_mle_rerun_removes_stale_optional_details(tmp_path, monkeypatch):
    args, _, _ = mle_fixture(tmp_path)
    def success(parsedargs):
        for suffix in ('.gene_summary.txt', '.sgrna_summary.txt'):
            Path(parsedargs.output_prefix + suffix).write_text('new ' + suffix)
        return ({}, None)
    monkeypatch.setattr('mageck2.mlemageck.mageckmle_main', success)
    run(args)
    assert not (tmp_path / 'result.fit-details.json').exists()
    manifest = json.loads((tmp_path / 'result.crisprworks.json').read_text())
    assert manifest['status'] == 'complete'
    assert set(manifest['outputs']) == {'gene_summary', 'sgrna_summary'}
    assert manifest['options']['output_prefix'] == str(tmp_path / 'result')


def test_mle_changed_input_aborts_publication(tmp_path, monkeypatch):
    args, counts, previous = mle_fixture(tmp_path)
    def changing_fit(parsedargs):
        assert Path(parsedargs.count_table).read_bytes() == counts.read_bytes()
        counts.write_text('changed during fitting')
        for suffix in ('.gene_summary.txt', '.sgrna_summary.txt'):
            Path(parsedargs.output_prefix + suffix).write_text('new result')
        return ({}, None)
    monkeypatch.setattr('mageck2.mlemageck.mageckmle_main', changing_fit)
    with pytest.raises(ValueError, match='Input changed'):
        run(args)
    assert all(path.read_bytes() == content for path, content in previous.items())


def test_failure_between_mle_table_replacements_leaves_no_completion_manifest(tmp_path, monkeypatch):
    files = {}
    for name in ('genes', 'guides'):
        output, staged = tmp_path / name, tmp_path / (name + '.tmp')
        output.write_text('old')
        staged.write_text('new')
        files[name] = (staged, output)
    manifest = tmp_path / 'manifest'
    manifest.write_text('old completion')
    replace = Path.replace
    def fail_guides(path, target):
        if target == tmp_path / 'guides':
            raise OSError('injected second replacement failure')
        return replace(path, target)
    monkeypatch.setattr(Path, 'replace', fail_guides)
    with pytest.raises(OSError, match='injected'):
        publish_file_bundle(files, manifest, {'status': 'complete'})
    assert not manifest.exists()
    assert (tmp_path / 'genes').read_text() == 'new'
    assert (tmp_path / 'guides').read_text() == 'old'


def test_mle_snapshots_preserve_csv_parsing(tmp_path, monkeypatch):
    from mageck2.mleinstanceio import read_gene_from_file
    args, counts, _ = mle_fixture(tmp_path)
    csv_counts = tmp_path / 'counts.csv'
    csv_counts.write_text(counts.read_text().replace('\t', ','))
    args.count_table = str(csv_counts)
    def parsed_fit(parsedargs):
        assert Path(parsedargs.count_table).suffix == '.csv'
        genes = read_gene_from_file(parsedargs.count_table)
        assert float(genes['A'].nb_count[0, 0]) == 11  # upstream adds one pseudocount
        for suffix in ('.gene_summary.txt', '.sgrna_summary.txt'):
            Path(parsedargs.output_prefix + suffix).write_text('new result')
        return ({}, None)
    monkeypatch.setattr('mageck2.mlemageck.mageckmle_main', parsed_fit)
    run(args)
    assert args.count_table == str(csv_counts)
    assert args.design_matrix == '1,0;1,1'


def test_optional_details_cleanup_cannot_delete_a_count_input(tmp_path):
    args, counts, _ = mle_fixture(tmp_path)
    details = tmp_path / 'result.fit-details.json'
    original = counts.read_bytes()
    details.write_bytes(original)
    args.count_table = str(details)
    assert not args.write_fit_details
    with pytest.raises(ValueError, match='collides'):
        run(args)
    assert details.read_bytes() == original


@pytest.mark.parametrize('value', ['nan', 'inf', '-1', 'not-a-count'])
def test_invalid_mle_counts_fail_before_upstream_and_preserve_prior_bundle(tmp_path, monkeypatch, value):
    args, counts, previous = mle_fixture(tmp_path)
    counts.write_text(f'sgRNA\tGene\tbaseline\ttreated\na\tA\t{value}\t5\n')
    def unexpected_fit(**kwargs):
        pytest.fail('Malformed counts reached the upstream fitter')
    monkeypatch.setattr('mageck2.mlemageck.mageckmle_main', unexpected_fit)
    with pytest.raises(ValueError, match='finite and nonnegative'):
        run(args)
    assert all(path.read_bytes() == content for path, content in previous.items())


def test_duplicate_mle_guides_fail_before_upstream(tmp_path, monkeypatch):
    args, counts, previous = mle_fixture(tmp_path)
    with counts.open('a') as stream:
        stream.write('a\tOTHER_GENE\t4\t5\n')
    monkeypatch.setattr('mageck2.mlemageck.mageckmle_main', lambda **kwargs: pytest.fail('Duplicate guide reached fitting'))
    with pytest.raises(ValueError, match='Duplicate guide'):
        run(args)
    assert all(path.read_bytes() == content for path, content in previous.items())
