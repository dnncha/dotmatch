# Temporary EditWitness delivery branch

This is an isolated delivery branch in dnncha/dotmatch, not part of DotMatch. Do not merge it into main. The project belongs in its own repository, dnncha/editwitness, which has not been created by the current connection.

The checked source transport holds 68 UTF-8 source files, including a SHA-256 inventory. It contains code, tests, synthetic examples, schemas, scientific documentation, a structured roadmap, human/agent instructions, license, and an independent-repository publishing utility. The bootstrap workflow is restricted to this exact branch and repository. It never targets main and does not publish to PyPI.

After the bootstrap materializes the source:

```bash
python scripts/release_manifest.py --check
python scripts/publish_github.py --public --dry-run
python scripts/publish_github.py --public
```

The last command requires locally installed and authenticated `gh` and `git`. It creates a fresh source-only repository, excluding this transport, temporary workflow, and DotMatch history. It refuses to overwrite an existing repository. Review BUILD_STATUS.md and the actual workflow run before asserting any validation gate passed.

Compressed source transport SHA-256:
`8ed72d20468f41f48dabcd1b6b93ee89e20076195e9800ba21ff6d76767a3486`

This integrity checksum detects accidental transport changes; it is not an authenticated signature. EditWitness is a research alpha, not an empirically validated assay or clone-safety certificate.
