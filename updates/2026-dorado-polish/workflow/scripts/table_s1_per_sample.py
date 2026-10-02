"""Table S1: per-Sample results for every Read set and Arm (#14).

One row per Sample x Read model x Depth x Arm, with the Read set's actual depth (overall and
per contig, from depth.tsv), then for SNP and INDEL the Best F1 (the QUAL sweep's best
threshold, with the QUAL threshold it was reached at) and the Default-PASS score, each with its
precision and recall. A CSV with readable headers, for the site's interactive csv-table.
Numbers are passed through as results.tsv wrote them.
"""

import sys

import pandas as pd

sys.stderr = open(snakemake.log[0], "w")

labels = snakemake.params.arm_labels
results = pd.read_csv(snakemake.input.results, sep="\t", dtype=str)
depth = pd.read_csv(snakemake.input.depth, sep="\t", dtype={"mean_depth": float})

KEYS = ["sample", "read_model", "depth"]
# "chromosome 5.03; plasmid 5.13" for each Read set, in the order the contigs are in depth.tsv
by_contig = (
    depth.assign(text=lambda d: d["contig"] + " " + d["mean_depth"].map("{:.2f}".format))
    .groupby(KEYS, sort=False)["text"]
    .agg("; ".join)
    .rename("actual_depth_by_contig")
    .reset_index()
)
by_contig["depth"] = by_contig["depth"].astype(str)

METRICS = [("precision", "precision"), ("recall", "recall"), ("f1", "F1")]
wide = results[results["var_type"].isin(["SNP", "INDEL"])]
table = (
    wide[["sample", "species", "dnd_sample", "read_model", "depth", "actual_depth", "arm"]]
    .drop_duplicates()
    .merge(by_contig, on=KEYS, how="left", validate="many_to_one")
)
columns = {}
for var_type in ["SNP", "INDEL"]:
    for mode, mode_name in [("sweep_best", "Best"), ("default_pass", "Default-PASS")]:
        part = wide[(wide["var_type"] == var_type) & (wide["scoring_mode"] == mode)]
        part = part.set_index(KEYS + ["arm"])
        for column, name in METRICS:
            columns[f"{var_type} {mode_name} {name}"] = (part, column)
        if mode == "sweep_best":
            columns[f"{var_type} Best QUAL threshold"] = (part, "min_qual")
for name, (part, column) in columns.items():
    table[name] = part[column].reindex(pd.MultiIndex.from_frame(table[KEYS + ["arm"]])).to_numpy()

table["arm_label"] = table["arm"].map(lambda a: f"{a} ({labels[a]})")
table = table.rename(
    columns={
        "sample": "Sample",
        "species": "Species",
        "dnd_sample": "dnd Sample",
        "read_model": "Read model",
        "depth": "Depth (x)",
        "actual_depth": "Actual depth (x)",
        "actual_depth_by_contig": "Actual depth by contig (x)",
        "arm_label": "Arm",
    }
)
table["Depth (x)"] = table["Depth (x)"].astype(int)
table = table.sort_values(["Sample", "Read model", "Depth (x)", "arm"], kind="stable")
first = ["Sample", "Species", "dnd Sample", "Read model", "Depth (x)", "Actual depth (x)",
         "Actual depth by contig (x)", "Arm"]  # fmt: skip
table = table[first + [c for c in table.columns if c not in first and c != "arm"]]
table.to_csv(snakemake.output.csv, index=False)
print(f"rows={len(table)} columns={len(table.columns)}", file=sys.stderr)
