import re
import sys
from pathlib import Path

sys.stderr = open(snakemake.log[0], "w")


field_names = [
    ("Number of alignments", "num_alignments", int),
    ("% from total reads", "percent_total_reads", float),
    ("Yield [Gb]", "yield_gb", float),
    ("Mean coverage", "mean_coverage", float),
    ("Yield [Gb] (>25kb)", "yield_gb_gt_25kb", float),
    ("N50", "n50", int),
    ("N75", "n75", int),
    ("Median length", "median_length", float),
    ("Mean length", "mean_length", float),
    ("Median identity", "median_identity", float),
    ("Mean identity", "mean_identity", float),
]
value_re = re.compile(r"(\d+(\.\d+)?)$")


def parse_stats_file(stats_file: Path):
    """Parse the stats file, handling the fact that they could appear anywhere in the file, and return a dict of the values"""
    stats = {}
    # name is filename without suffix
    name = stats_file.name.rsplit(".", maxsplit=1)[0]
    model = stats_file.parent.name
    stats["filename"] = name
    stats["model"] = model
    with open(stats_file) as fh:
        for line in map(str.strip, fh):
            for field_name, new_field_name, dtype in field_names:
                if field_name in line:
                    stats[new_field_name] = dtype(value_re.search(line).group(1))
    return stats


def main():
    out_fp = open(snakemake.output.stats, "w")

    stats = {}
    for stats_file in map(Path, snakemake.input.stats):
        stats[stats_file] = parse_stats_file(stats_file)

    header_written = False
    for p, d in sorted(stats.items()):
        if header_written is False:
            header_written = True
            print(",".join(d.keys()), file=out_fp)
        print(",".join(map(str, d.values())), file=out_fp)

    return 0


if __name__ == "__main__":
    sys.exit(main())
