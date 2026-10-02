"""Table 1: runtime and memory of each alignment and calling step (#14).

Built from benchmarks.tsv. A row is a step (alignment or variant calling) of one tool at one
Depth, aggregated over the Samples and Read models in the run: the median and range of wall
time (seconds) and peak memory (MB). Two things keep the rows honest:
  - Arms with the same aligner, version and preset share one alignment, so benchmarks.tsv
    repeats its row for each of them. Alignment rows are counted once per Read set and list
    all the Arms that use them.
  - Dorado's `--device cpu` re-run (timing only, never scored) is its own row, beside the
    GPU one, at the Depths in dorado.cpu_run.
Peak memory is the whole job's peak RSS, which is host memory: GPU memory isn't measured.
The headers are readable because the CSV is shown on the site as a csv-table.
"""

import sys

import pandas as pd

sys.stderr = open(snakemake.log[0], "w")

labels = snakemake.params.arm_labels
bench = pd.read_csv(snakemake.input.benchmarks, sep="\t", dtype={"tool_version": str})
bench["timing_only"] = bench["timing_only"].astype(str).str.lower() == "true"

TOOLS = {"minimap2": "minimap2", "clair3": "Clair3", "dorado": "dorado polish"}
STEPS = {"align": "Alignment", "call": "Variant calling"}

# One alignment per Read set and alignment, however many Arms use it.
align = bench[bench["step"] == "align"].copy()
align["arms"] = align.groupby(["sample", "read_model", "depth", "alignment"])["arm"].transform(
    lambda arms: ", ".join(sorted(set(arms)))
)
align = align.drop_duplicates(["sample", "read_model", "depth", "alignment"])
call = bench[bench["step"] == "call"].copy()
call["arms"] = call["arm"]
rows = pd.concat([align, call])

# What distinguishes one table row from another, besides the Depth.
GROUP = ["step", "tool", "tool_version", "alignment", "arms", "device", "timing_only", "hardware", "threads"]
out = []
for keys, grp in rows.groupby(GROUP + ["depth"]):
    key = dict(zip(GROUP + ["depth"], keys))
    wall, rss = grp["wall_time_s"], grp["max_rss_mb"]
    preset = key["alignment"].split("-", 2)[2] if key["step"] == "align" else ""
    arms = key["arms"]
    out.append(
        {
            "Step": STEPS[key["step"]],
            "Arms": ", ".join(f"{a} ({labels[a]})" for a in arms.split(", ")),
            "Tool": TOOLS[key["tool"]],
            "Version": key["tool_version"],
            "Alignment": key["alignment"],
            "Device": key["device"].upper(),
            "Timing only": "yes" if key["timing_only"] else "no",
            "Hardware": key["hardware"],
            "Threads": key["threads"],
            "Depth (x)": key["depth"],
            "Read sets": len(grp),
            "Wall time median (s)": round(wall.median(), 2),
            "Wall time min (s)": round(wall.min(), 2),
            "Wall time max (s)": round(wall.max(), 2),
            "Peak RSS median (MB)": round(rss.median(), 1),
            "Peak RSS min (MB)": round(rss.min(), 1),
            "Peak RSS max (MB)": round(rss.max(), 1),
            "_order": (list(STEPS).index(key["step"]), arms, key["timing_only"], key["depth"]),
        }
    )
table = pd.DataFrame(out)
table = table.sort_values("_order").drop(columns="_order")
table.to_csv(snakemake.output.csv, index=False)
print(f"rows={len(table)} read_sets_per_row={sorted(set(table['Read sets']))}", file=sys.stderr)
