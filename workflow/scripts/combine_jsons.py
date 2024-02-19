"""This script combines a long list of JSONs into a single JSON."""

import sys

sys.stderr = open(snakemake.log[0], "w")
import json


def main():
    jsons = []
    for json_file in snakemake.input.jsons:
        with open(json_file) as f:
            jsons.append(json.load(f))

    with open(snakemake.output.json, "w") as out:
        json.dump(jsons, out, indent=4)


if __name__ == "__main__":
    main()
