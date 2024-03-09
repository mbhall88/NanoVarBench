"""take a fastq file (from stdin) and a pod5 read summary file and output the read ids in the fastq 
that exist in the fastq file, accounting for the fact that due to read splitting,
some fastq read ids don't exist in the pod5s and thus have a pi:Z tag with the correct id
"""

import sys
import pysam
import re

PI_REGEX = re.compile(r"\spi:Z:(?P<oid>[A-Za-z0-9-]+)\s")


def main():
    pod5_summary = sys.argv[1]
    with open(pod5_summary) as f:
        _ = next(f)
        pod5_ids = set(line.split(",")[0].strip() for line in f)

    print(f"Read {len(pod5_ids)} read ids from {pod5_summary}", file=sys.stderr)

    output_ids = set()

    for read in pysam.FastxFile("-"):
        read_id = read.name
        pi_match = PI_REGEX.search(read.comment)
        if pi_match:
            original_id = pi_match.group("oid")
            if original_id in pod5_ids:
                if read_id in pod5_ids:
                    raise ValueError(
                        f"A split read's original ID and subread ID are both in the pod5s: subread: {read_id} and original: {original_id}"
                    )
                else:
                    if original_id not in output_ids:
                        print(original_id)
                        output_ids.add(original_id)
            else:
                if read_id not in output_ids:
                    raise ValueError(
                        f"Original ID {original_id} and subread ID {read_id} not in pod5s"
                    )
        else:
            if read_id in pod5_ids:
                if read_id not in output_ids:
                    print(read_id)
                    output_ids.add(read_id)
            else:
                raise ValueError(f"Read ID {read_id} not in pod5s")


if __name__ == "__main__":
    main()
