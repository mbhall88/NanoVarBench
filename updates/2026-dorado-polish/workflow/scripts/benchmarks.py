"""Build the benchmark table (Table 1, #12) from Snakemake's benchmark files.

benchmarks.tsv has one row per Sample x Read model x Depth x Arm x step x timing_only:
  step         align (the Arm's alignment, so Arms sharing a BAM report the same alignment
               cost) or call (Clair3 or `dorado polish`)
  timing_only  true for the Dorado CPU re-run at the Depths in dorado.cpu_run, a second call
               row for Arm D. These runs are never scored: they aren't in results.tsv.
  device       gpu or cpu: where the rule runs, from the config.
  hardware     the config's benchmark_hardware for that device (the GPU model, or
               "CPU: <model>"). Nothing is recorded in the jobs: the Slurm profile pins the
               partition and constraint that provide this hardware.
  threads      the rule's configured threads.
wall_time_s is Snakemake's `s`, the whole job script, which only runs the tool (plus sorting
and indexing for align). Memory is in MB and times in seconds, as in Snakemake's benchmark
files.
"""

import csv
import sys

sys.stderr = open(snakemake.log[0], "w")

# Snakemake benchmark column -> our column. max_rss, max_vms, max_uss, max_pss, io_in and io_out
# are in MB, s and cpu_time in seconds, mean_load is a CPU percentage.
BENCHMARK_COLUMNS = {
    "s": "wall_time_s",
    "max_rss": "max_rss_mb",
    "max_vms": "max_vms_mb",
    "max_uss": "max_uss_mb",
    "max_pss": "max_pss_mb",
    "cpu_time": "cpu_time_s",
    "mean_load": "mean_cpu_pct",
    "io_in": "io_in_mb",
    "io_out": "io_out_mb",
}
COLUMNS = [
    "sample", "read_model", "depth", "arm", "step", "timing_only", "tool", "tool_version",
    "alignment", "device", "hardware", "threads",
    "wall_time_s", "max_rss_mb", "max_vms_mb", "max_uss_mb",
    "max_pss_mb", "cpu_time_s", "mean_cpu_pct", "io_in_mb", "io_out_mb",
]  # fmt: skip
FROM_JOB = [
    "sample", "read_model", "depth", "arm", "step", "timing_only", "tool", "tool_version",
    "alignment", "device", "threads",
]  # fmt: skip


def read_one(path):
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if len(rows) != 1:
        raise ValueError(f"{path}: expected one row, found {len(rows)}")
    return rows[0]


rows = []
for job in snakemake.params.jobs:
    bench = read_one(job["benchmark"])
    rows.append(
        {
            **{k: job[k] for k in FROM_JOB},
            "hardware": snakemake.params.hardware[job["device"]],
            **{ours: bench[theirs] for theirs, ours in BENCHMARK_COLUMNS.items()},
        }
    )

with open(snakemake.output.tsv, "w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
print(f"benchmarks={len(rows)}", file=sys.stderr)
