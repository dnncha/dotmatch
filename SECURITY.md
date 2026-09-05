# Security and sensitive data

EditWitness's scientific analysis runs locally and performs no network access,
telemetry, reference fetching or arbitrary plugin loading. Installing dependencies
can use the configured package index; use your institution's approved environment
and dependency controls.

Full JSON and HTML may contain genomic sequence. Store and share them deliberately.
A self-contained report is not anonymized merely because it needs no server.
Error messages avoid echoing rejected Pydantic input values, but paths and user
identifiers can appear. Do not put sensitive identifiers in public bug reports.

HTML escapes supplied text, uses no JavaScript or remote assets, and has a
restrictive content-security policy. Input sizes, sequence work, scan grids and
HTML previews are bounded. Files are created without silent replacement and
are individually atomic; a multi-output command is not an atomic transaction.

A SHA-256 result checksum detects accidental changes when the expected checksum
is trusted. It is not a cryptographic author signature, a secure execution
attestation, or biological validation. Anyone can recompute a checksum.

Until a standalone private vulnerability-reporting channel is configured,
contact the maintainer through their existing private GitHub contact route rather
than publishing exploitable details or private sequences. The release checklist
requires enabling GitHub private vulnerability reporting where available.
