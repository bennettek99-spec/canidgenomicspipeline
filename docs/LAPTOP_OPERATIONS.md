# Laptop operation and recovery

CANIS defaults to a conservative local operating model. It is not a cluster scheduler: the intended workflow is reduced panels or modest pre-called VCFs, sequential chromosome work, bounded workers, and resumable validated stages.

## Resource manager

The `resource_manager` configuration block is enforced by the native executor.

```yaml
executor:
  max_workers: 2
resource_manager:
  max_workers: 2
  max_memory_mb: 12288
  min_free_disk_mb: 8192
  temperature_friendly: true
  max_threads_per_stage: 4
  process_chromosomes_sequentially: true
  cleanup_intermediates: true
```

Before a run, CANIS checks available disk space across the data, run, and cache locations. It limits independent stage batches by declared memory, leaves thread headroom for a thin laptop, sets conservative numerical-library thread defaults if they are not already set, and removes stale temporary transactions only when configured to do so.

## Pause, stop, and resume

`canidae ui` exposes buttons for safe controls. The command-line equivalent is to create markers in the chosen run directory:

- `runs/<run-id>/.canidae.pause` pauses at the next stage-batch boundary.
- `runs/<run-id>/.canidae.cancel` stops at the next stage-batch boundary.

Remove the marker and rerun with the same `--run-id` to continue. The safe cache accepts a previous stage only when its inputs, config, code fingerprint, and registered output hashes still agree. The run manifest includes the recovery command context and a reason for every skip or stop.

## Output lifecycle

Each stage writes under a temporary `.staging/` directory. After all declared outputs exist and can be hashed, CANIS promotes that stage directory atomically and atomically updates `artifacts.json`. A crash therefore leaves either the previous completed result or no registered new result; it cannot expose a half-written VCF as a completed artifact.

The cache contains JSON fingerprints only, never duplicate genomics files. Use `executor.resume: false` to force recalculation, or remove an individual cache record under `data/cache/stage-cache/` when intentionally rebuilding one stage.

## Suggested laptop recipe

1. Begin with the 2K preset to test connectivity and metadata.
2. Use 10K for a normal exploratory analysis; 25K only after reviewing the real transfer preflight.
3. Keep `max_workers` at 1–2 for a thin-and-light machine; let the UI or laptop profile set memory and disk guardrails.
4. Open `data/store/report/report.html` directly after a successful run. Treat its panel-relative warnings and methods section as part of the result, not footnotes.

Cross-assembly liftover, high-coverage raw-read alignment, and full callable-genome demography are deliberate preparation or larger-compute workflows, not automatic laptop jobs.
