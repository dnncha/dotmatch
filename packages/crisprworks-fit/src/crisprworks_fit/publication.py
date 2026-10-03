"""Exclusive staging files beside their final destination."""

from pathlib import Path
import os
import tempfile
import hashlib
import json
from contextlib import contextmanager


def file_digest(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_text_input(path):
    """Parseable text and provenance from exactly the same immutable bytes."""
    path = Path(path)
    payload = path.read_bytes()
    return payload.decode("utf-8"), {
        "path": str(path.resolve()), "sha256": hashlib.sha256(payload).hexdigest(),
    }


def verify_output_record(metadata):
    """Verify a recorded companion table; legacy standalone details have none."""
    if "output" not in metadata:
        return
    record = metadata["output"]
    if (not isinstance(record, dict) or not isinstance(record.get("path"), str)
            or not record["path"] or not isinstance(record.get("sha256"), str)
            or len(record["sha256"]) != 64
            or any(c not in "0123456789abcdef" for c in record["sha256"])):
        raise ValueError("Malformed companion output record")
    if file_digest(record["path"]) != record["sha256"]:
        raise ValueError("Companion output hash mismatch; fit bundle is inconsistent")


def reject_input_collisions(destinations, inputs):
    """Protect path aliases and distinct hard links to the same input inode."""
    inputs = [Path(value) for value in inputs if value is not None]
    for destination in map(Path, destinations):
        for source in inputs:
            if destination.resolve() == source.resolve() or (
                destination.exists() and source.exists() and destination.samefile(source)
            ):
                raise ValueError("Output collides with an input file")


def unique_staging_path(destination):
    """Reserve a new file; never reuse a predictable input or another run's file."""
    destination = Path(destination)
    descriptor, name = tempfile.mkstemp(
        prefix=f".{destination.name}.", suffix=".tmp", dir=destination.parent,
    )
    os.close(descriptor)
    return Path(name)


@contextmanager
def publication_lock(manifest):
    """Fail on overlapping publication; the OS releases the lock after a crash.

    Keep the lock file in place: unlinking it would let another process lock
    a different inode while the first publisher still holds its descriptor.
    """
    path = manifest.with_name(f".{manifest.name}.publish.lock")
    descriptor = os.open(path, os.O_CREAT | os.O_RDWR | getattr(os, "O_NOFOLLOW", 0), 0o600)
    try:
        if os.name == "nt":
            import msvcrt
            if os.fstat(descriptor).st_size == 0:
                os.write(descriptor, b"\0")
            os.lseek(descriptor, 0, os.SEEK_SET)
            try:
                msvcrt.locking(descriptor, msvcrt.LK_NBLCK, 1)
            except OSError as error:
                raise ValueError(f"Publication already in progress for {manifest}") from error
        else:
            import fcntl
            try:
                fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as error:
                raise ValueError(f"Publication already in progress for {manifest}") from error
        yield
    finally:
        os.close(descriptor)


def verify_input_records(metadata):
    for name, record in metadata.get("inputs", {}).items():
        if isinstance(record, dict) and "path" in record:
            if file_digest(record["path"]) != record.get("sha256"):
                raise ValueError(f"Input changed before publication: {name}")


def publish_table_bundle(staged_table, table, manifest, metadata):
    """Stage everything before replacement; manifest absence marks an incomplete run.

    This protects against staging failures and stale manifests on interrupted
    publication. Overlapping cooperating publishers fail before replacement.
    """
    metadata["output"] = {
        "path": str(table.parent.resolve() / table.name), "sha256": file_digest(staged_table),
    }
    _publish_bundle([(staged_table, table)], manifest, metadata)


def publish_file_bundle(outputs, manifest, metadata, remove=()):
    """Publish all required MLE files with the completion manifest last."""
    metadata["outputs"] = {
        name: {"path": str(destination.parent.resolve() / destination.name),
               "sha256": file_digest(staged)}
        for name, (staged, destination) in outputs.items()
    }
    _publish_bundle(list(outputs.values()), manifest, metadata, remove)


def _publish_bundle(files, manifest, metadata, remove=()):
    staged_manifest = unique_staging_path(manifest)
    destinations = [destination for _, destination in files] + list(remove) + [manifest]
    try:
        staged_manifest.write_text(json.dumps(metadata, indent=2, allow_nan=False) + "\n")
        inputs = [record["path"] for record in metadata.get("inputs", {}).values()
                  if isinstance(record, dict) and "path" in record]
        lock_path = manifest.with_name(f".{manifest.name}.publish.lock")
        reject_input_collisions((*destinations, lock_path), inputs)
        with publication_lock(manifest):
            # Recheck after staging, immediately before replacing destinations.
            reject_input_collisions(destinations, inputs)
            verify_input_records(metadata)
            # Do not leave an earlier completion record beside a new table.
            manifest.unlink(missing_ok=True)
            for destination in remove:
                destination.unlink(missing_ok=True)
            for staged, destination in files:
                staged.replace(destination)
            staged_manifest.replace(manifest)
    finally:
        staged_manifest.unlink(missing_ok=True)
