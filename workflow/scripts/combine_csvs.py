"""This script combines a long list of CSVs into a single CSV, keeping only the header 
from the first file.
Assumptions:
- All CSVs have the same header
- All CSVs have the same number of columns
- All CSVs have the same delimiter
"""

import sys

sys.stderr = open(snakemake.log[0], "w")
import fileinput


def main():
    with open(snakemake.output.csv, "w") as out:
        input_csvs = snakemake.input.csvs
        header_written = False
        for line in fileinput.input(input_csvs):
            if fileinput.isfirstline():
                if not header_written:
                    out.write(line)
                    header_written = True
                continue
            out.write(line)


if __name__ == "__main__":
    main()
