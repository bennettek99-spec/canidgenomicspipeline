# Laptop-safe FASTQ archive acquisition

`acquire_reads` is the guarded raw-read counterpart to the indexed `reduced_panel` VCF extractor. It supports public gzip FASTQ archives in three layouts:

| Layout | URLs per sample | Downsampling behavior | Alignment handoff |
| --- | ---: | --- | --- |
| `paired` | 2 (R1/R2) | One random decision per mate pair, so mates never diverge | `bwa-mem2 mem R1 R2` |
| `interleaved` | 1 | One random decision per adjacent two-record pair | `bwa-mem2 mem -p` |
| `single` | 1 | One random decision per FASTQ record | `bwa-mem2 mem` |

## Safety contract

Before requesting archive bodies, CANIS reads `expected_bytes` from configuration or queries HTTP `Content-Length`/a one-byte range response. It sums the estimate across all selected archives and:

- refuses more than the hard 9 GB default ceiling;
- asks for `confirm_large_transfer: true` above the configured confirmation threshold;
- checks free disk space for the estimated retained download plus a safety overhead;
- streams compressed archive bytes through an MD5/SHA-256 verifier while decompressing FASTQ records;
- writes only `.part` files until an archive completes and its publisher checksum passes;
- removes `.part` files after every failure.

Publisher checksums should use `md5:<digest>` or `sha256:<digest>`. If an archive has no publisher checksum, CANIS still writes source SHA-256 and MD5 values to `read_acquisition_manifest.json`, but that is a recorded fingerprint rather than third-party verification.

## Recipe

Copy [the template](../configs/examples/laptop_read_acquisition.yaml), replace every `REPLACE_*` sample and URL, and update the sample sheet to include precisely those sample IDs. Use a small `downsample_fraction` first.

```powershell
canidae run -c configs/examples/laptop_read_acquisition.yaml `
  -c configs/profiles/laptop.yaml
```

`acquire_reads` writes a `fastq_manifest.csv` artifact with `sample_id`, `layout`, `fastq1`, and `fastq2`. The normal `align` stage detects that artifact automatically; no manual file renaming or paired/interleaved conversion is required.

Raw-read calling remains a deliberate WSL/Linux workflow because `bwa-mem2`, `samtools`, and `bcftools` are external bioinformatics tools. For ordinary laptop exploration, prefer the compact indexed VCF [reduced-panel workflow](REDUCED_PANEL_PIPELINE.md).
