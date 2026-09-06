# Migration from 0.1.0a1

Version 0.2.0a1 emits result schema 1.1. It accepts manifest schemas 1.0 and 1.1.

**An omitted `observation_model` still means `original-sites-presence-v1`.** To opt into sequence-aware rematching, add:

```json
{"schema_version": "1.1", "observation_model": "exact-local-sites-presence-v2"}
```

These keys are additions to a full manifest, not a complete manifest by themselves. `demo` and `init` now make this selection explicitly. `demo --legacy-model` uses the historical behavior. `compare-models` shows whether witnesses or selected panels change.

Analysis results now contain the input allele edit definitions, the number of different local genotypes, and each hypothesis's representative. Sequence-equivalent hypotheses share one witness. Two new conclusions prevent empty comparisons or negative baselines from looking reassuring: `no_distinct_alternatives` and `baseline_uninformative`.

Exact observations can contain several `products`. Each product includes its edited-sequence coordinates, F-primer-oriented readout and signal ID. The old singular `reads`, `product_length` and `signal_id` fields are populated only for a single product; they are not the authoritative signal set for a multiple-product result. Use `hypothesis_observations.signal_ids` for genotype comparisons or the products collection for detailed evidence. Product coordinates and latent product lengths are diagnostic metadata, not extra measured features.

The `scan` command remains a geometry-only original-site calculation. Its model identifier does not change to the manifest's sequence-aware selection.

Use the producing 0.1.0a1 package to verify old result-schema 1.0 evidence. New code explicitly rejects unsupported result schemas; it does not rewrite and rehash old scientific results as though they were unchanged.
