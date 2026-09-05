# Roadmap

The initial engine, CLI, schemas, local report, examples, tests and packaging are
implemented. The next work must increase scientific usefulness and adoption,
not expand the number of unsupported feature claims.

The authoritative task list is [roadmap.json](../roadmap.json). Tasks include
acceptance criteria, dependencies, scope and status so a later agent can continue
without another architecture exercise.

## Priority order

**First: independent scrutiny and reliable publication.** Finish the standalone
repository/package-index setup, execute the full CI matrix, and have an
independent genome-engineering scientist challenge the observation model. Publish
as a research alpha, not as experimentally validated software.

**Next: sequence-aware primer-site modeling.** The original-site model is explicit
but can be overly conservative or representation-dependent for replacements
that recreate binding sequence. Add a separately versioned exact edited-sequence
model, ambiguous-product handling, and adversarial normalization tests before
adding thermodynamic mismatch probabilities.

**Then: independently adjudicated examples and workflow adapters.** Collect real
reference/primer/ground-truth metadata with permission and provenance. Add a
read-only adapter for one established caller format only after its semantics
are pinned and tested. An adapter must never infer that an unobserved outcome has
zero biological frequency.

**Then: useful orthogonal measurements.** A copy-number-aware response model needs
measurement uncertainty and calibration, not an idealized integer pasted onto
the current observation set. Keep those claims gated until that evidence exists.

**Finally: adoption-driven breadth.** More hypothesis generators, batch workflow
examples, Conda/Bioconda packaging, workflow registries and a thin optional MCP
adapter should follow independent use. Profile before introducing Rust, GPUs or
distributed computation.

## Work we will not quietly do

No invented public benchmark results. No clinical certification. No confidence
scores without a justified statistical model. No "AI-powered" branding for a
deterministic engine. No opaque auto-generated hypothesis set presented as
biological completeness. No replacement of trusted aligners merely to increase
this project's scope.
