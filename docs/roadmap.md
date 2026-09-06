# Roadmap

`../roadmap.json` is the machine-readable source of task status and acceptance
criteria. A status of implemented is not empirical biological validation.

The 0.2 audit delivered separately versioned exact primer rematching, bounded
reference-deletion hypothesis generation, canonical genotype grouping and model
comparison. Those capabilities now need independent scientific review and real
cases, rather than more speculative features.

The next highest-value scientific work is EW-002 and EW-004: external review of
both observation models and an adjudicated public benchmark with exact assay
metadata and independently established outcomes. No reviewer or dataset result
has been fabricated. Include negative controls and false assurances, not only
cases where the program produces a useful warning.

The next distribution work is EW-001: move the public independent source into
its own repository, establish package-index ownership, then configure reviewed
trusted publishing. Public GitHub source, a GitHub prerelease and PyPI publication
are different states and must be reported separately.

Only then choose a real read-only caller adapter, calibrate an orthogonal
measurement model, and evaluate repeat use at unrelated facilities. Preserve a
local, deterministic engine with transparent limitations throughout.
