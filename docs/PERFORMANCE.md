# Laptop performance and resource policy

CANIS is optimized for laptop-safe exploratory VCF analysis. The built-in default is not a large-cluster scheduler and it does not automatically launch raw-read downloads, liftover, or full-genome workflows.

| Workload | Expected laptop behavior |
| --- | --- |
| Indexed reduced panel (2K) | Small connectivity/QC smoke test; minutes after range preflight |
| Indexed reduced panel (10K) | Normal exploratory PCA/distance/tree workflow; tens of minutes |
| Indexed reduced panel (25K) | Higher-resolution exploratory workflow; review exact transfer estimate first |
| Same-build local VCF | Depends on VCF size; QC/load once, then compact genotype matrix reused |
| Local ancestry / D/f statistics | Chromosome-sequential to bound peak memory; slightly longer but stable |
| Raw reads / full calling | External Linux/WSL preparation workflow, not a default laptop command |

## Recommended profile

```powershell
canidae run -c configs/examples/redwolf_jackal_reduced_panel.yaml `
  -c configs/profiles/laptop.yaml
```

The laptop profile limits independent stages to two, reserves a memory ceiling, requires an 8 GiB free-disk floor, keeps numerical-library thread defaults at four or fewer, cleans stale transaction workspaces, and processes chromosome-aware analyses sequentially.

Adjust only after checking your hardware:

```yaml
executor:
  max_workers: 2
resource_manager:
  max_workers: 2
  max_memory_mb: 12288
  min_free_disk_mb: 8192
  max_threads_per_stage: 4
```

Use `canidae ui` for a browser view of these settings, preliminary reduced-panel estimates, and safe controls. The remote extractor computes the real compressed byte-range estimate before source VCF ranges are downloaded; its default hard ceiling is 9 GB.

See [Laptop operation and recovery](LAPTOP_OPERATIONS.md) for cache, pause/resume, cancellation, and cleanup behavior.
