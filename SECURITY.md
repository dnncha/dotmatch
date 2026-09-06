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

A dedicated private vulnerability-reporting channel is not yet configured.
Do not publish exploitable details, confidential sequences or personal data in
public issues, and do not send them until a private reporting route is confirmed.
Enabling GitHub private vulnerability reporting is a post-publication task.
Synthetic, non-sensitive correctness reproducers can be filed publicly.
