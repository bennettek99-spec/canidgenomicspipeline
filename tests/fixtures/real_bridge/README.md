# Real bridge-panel regression fixture

This is a committed, derived subset of the public bridge-panel inputs used by
the NYC coydog and eastern-coyote recipes. It contains 252 bridge loci, four
query call tables, the corresponding reference calls, and 31 WGS reference
samples spanning coyote, wolf, village-dog, and breed groups.

Sources are the public inputs prepared under `data/`:

- `eastern_coyote_wgs_bridge/eastern_coyote_wgs_bridge_1000.vcf.gz`
- `nyc_coydog_validation/`
- `nyc_coydog_breeds/all_sample_genotypes.json`

The full study panels remain outside Git and are covered by the local promotion
tests. This fixture exists only to keep a small real-data regression surface in
CI; it is not a replacement for the full panels or a new scientific dataset.
