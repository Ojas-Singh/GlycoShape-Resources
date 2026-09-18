# ReGlyco Colab workflows

The two notebooks in this directory run the native Rust `reglyco` executable
locally in the Colab VM. GlycoShape remains the level-1 asset provider used by
the executable; the provider setting is fixed in `reglyco_local.py` and is not
an interactive notebook option.

`reglyco_local.py` resolves a pinned Linux x86-64 release described by
`reglyco_release.json`. Until that GitHub release asset is uploaded, it falls
back to a checkout selected with `REGLYCO_SOURCE_ROOT` or to:

```bash
cargo install reglyco --locked --version 0.1.0
```

For local testing, set `REGLYCO_BIN` to an executable path. The optional
`REGLYCO_BINARY_URL` and `REGLYCO_SHA256` variables override the release URL
and checksum. Every command writes stdout/stderr logs and records its command,
seed, executable checksum, provider, and output paths in the notebook run
artifacts.

The binder-design notebook uses local `scan`, `ensemble`, and `build` commands
and keeps its existing Pass/Borderline/Fail and Mol* inspection workflow. The
local-ensemble notebook treats its requested count as accepted glycosylated
structures, attaches `G00028MO` at the selected target residue, retains
rejection records, and writes a separate glycosylated multiframe PDB so the
glycan atoms remain visible in Mol*.
