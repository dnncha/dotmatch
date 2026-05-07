# Security Policy

DotMatch processes local sequence files and benchmark artifacts. It is not intended to run untrusted files as a network service without additional sandboxing.

## Supported Versions

Security fixes target the latest development branch until formal releases begin. After `v1.0`, this policy should be updated with supported release lines.

## Reporting a Vulnerability

Please do not open a public issue for a suspected vulnerability. Use GitHub private vulnerability reporting when it is enabled for the repository. If that is unavailable, contact the maintainer through the contact method listed on the GitHub profile or package metadata.

Useful reports include:

- a minimal input file or command that triggers the issue;
- operating system, compiler, and DotMatch commit;
- whether AddressSanitizer or UndefinedBehaviorSanitizer reports an error;
- whether the issue affects `dotmatch`, the Python bindings, Docker usage, or benchmark scripts.

## Scientific Data

Do not attach private sequencing data, patient data, human subject data, unpublished datasets, or proprietary sample sheets to public issues. Reduce reports to synthetic fixtures whenever possible.

