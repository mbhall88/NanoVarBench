import argparse
import gzip
import logging
import re
from collections import Counter
from enum import Enum
from pathlib import Path

import pysam

BARCODE_RE = re.compile(r"(barcode\d\d|unclassified)")


class ReadType(Enum):
    DUPLEX = 1
    SIMPLEX_WITH_OFFSPRING = -1
    SIMPLEX_NO_OFFSPRING = 0


def read_summary_file(summary_file: str) -> dict[str, str]:
    with open(summary_file) as fp:
        headers = next(fp).strip().split("\t")
        try:
            read_id_idx = headers.index("read_id")
        except ValueError:
            raise ValueError("read_id column not found in summary file")
        try:
            barcode_idx = headers.index("barcode")
        except ValueError:
            raise ValueError("barcode column not found in summary file")

        id2barcode = {}
        for line in fp:
            fields = line.strip().split("\t")
            read_id = fields[read_id_idx]
            barcode = fields[barcode_idx]
            id2barcode[read_id] = barcode

    return id2barcode


def extract_barcode(s: str) -> str:
    m = BARCODE_RE.search(s)
    if not m:
        raise ValueError(f"Invalid barcode: {s}")
    return m.group(1)


def to_fastq(read: pysam.AlignedSegment) -> str:
    qual = "".join(map(lambda x: chr(x + 33), read.query_qualities))
    seq = read.query_sequence
    tags = " ".join(f"{tag}={value}" for tag, value in read.tags)
    return f"@{read.query_name} {tags}\n{seq}\n+\n{qual}\n"


def main():
    parser = argparse.ArgumentParser(description="Duplex Read Demultiplexer")
    parser.add_argument("summary_file", help="Path to the sequencing summary file")
    parser.add_argument(
        "bam_file",
        help="Path to the [BS]AM file [default: stdin]. Stdin must be a BAM file.",
        default="-",
    )
    parser.add_argument(
        "-d", "--duplex", action="store_true", help="Keep only duplex reads"
    )
    parser.add_argument(
        "-o",
        "--outdir",
        help="Output barcodes files to this directory [default: demux]",
        default="demux",
        type=Path,
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Enable verbose logging"
    )
    parser.add_argument(
        "-U",
        "--greedy",
        action="store_true",
        help="If a duplex read's barcodes do not match, and one is unclassified, use the classified barcode",
    )
    parser.add_argument(
        "-f",
        "--fastq",
        action="store_true",
        help="Output fastq files instead of BAM files",
    )
    parser.add_argument(
        "-z",
        "--gzip",
        action="store_true",
        help="Compress output fastq files with gzip",
    )
    parser.add_argument(
        "-s",
        "--fastq-suffix",
        help="Suffix to add to fastq files [default: .fastq]",
        default=".fastq",
    )

    args = parser.parse_args()

    # Set logging level based on verbose flag
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG if args.verbose else logging.INFO)
    # Set logging format to include level name
    formatter = logging.Formatter(
        "[%(asctime)s] [%(levelname)s]: %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )
    handler = logging.StreamHandler()
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    summary_file = args.summary_file

    mode = "rb"
    if args.bam_file != "-" and args.bam_file.endswith(".sam"):
        mode = "r"
    # we don't check for sequences because these are generally unaligned BAMs
    bam_file = pysam.AlignmentFile(args.bam_file, mode=mode, check_sq=False)
    keep_duplex = args.duplex

    outdir = args.outdir.resolve()
    if outdir.exists() and not outdir.is_dir():
        raise ValueError(f"{outdir} exists and is not a directory")
    if not outdir.exists():
        outdir.mkdir(parents=True)

    handles = {}

    # Read the summary file
    id2barcode = read_summary_file(summary_file)
    logging.info(f"Read {len(id2barcode)} reads from {summary_file}")

    # Iterate over the BAM file and write the output
    num_mismatched_barcodes = 0
    num_orphan_reads = 0
    barcode_counts = Counter()

    logger.info("Demultiplexing reads...")
    for read in bam_file:
        try:
            dx = read.get_tag("dx")
        except KeyError:
            raise KeyError("dx tag not found in BAM file")

        read_type = ReadType(dx)

        if keep_duplex and read_type != ReadType.DUPLEX:
            continue

        if read_type != ReadType.DUPLEX:
            barcode = id2barcode.get(read.query_name)
            if barcode is None:
                logger.debug(
                    f"Barcode not found for simplex read: {read.query_name} - "
                    "this is a known issue https://github.com/nanoporetech/dorado/issues/474 - "
                    "marking as unclassified"
                )
                barcode = "unclassified"
                num_orphan_reads += 1
            barcode = extract_barcode(barcode)

        else:
            if ";" not in read.query_name:
                raise ValueError(
                    f"Duplex read IDs should have a ';' delimiter: {read.query_name}"
                )

            try:
                id1, id2 = read.query_name.split(";")
            except ValueError:
                raise ValueError(
                    f"Duplex read IDs should have exactly two read IDs: {read.query_name}"
                )

            barcode1 = id2barcode.get(id1)
            if barcode1 is None:
                logger.debug(
                    f"Barcode not found for duplex read: {id1} - "
                    "this is a known issue https://github.com/nanoporetech/dorado/issues/474 - "
                    "marking as unclassified"
                )
                barcode1 = "unclassified"
                num_orphan_reads += 1
            else:
                barcode1 = extract_barcode(barcode1)

            barcode2 = id2barcode.get(id2)
            if barcode2 is None:
                logger.debug(
                    f"Barcode not found for duplex read: {id2} - "
                    "this is a known issue https://github.com/nanoporetech/dorado/issues/474 - "
                    "marking as unclassified"
                )
                barcode2 = "unclassified"
                num_orphan_reads += 1

            barcode2 = extract_barcode(barcode2)

            if barcode1 != barcode2:
                if "unclassified" in (barcode1, barcode2) and args.greedy:
                    if barcode1 == "unclassified":
                        barcode = barcode2
                    else:
                        barcode = barcode1
                    logger.debug(
                        f"Barcodes do not match for duplex read: {read.query_name} - {barcode1} != {barcode2} - "
                        "but one is unclassified - using classified barcode"
                    )
                else:
                    logger.debug(
                        f"Barcodes do not match for duplex read: {read.query_name} - {barcode1} != {barcode2} - marking as unclassified"
                    )
                    num_mismatched_barcodes += 1
                    barcode = "unclassified"
            else:
                barcode = barcode1

        read.set_tag("BC", barcode, value_type="Z")

        fp = handles.get(barcode)
        if fp is None:
            path = str(outdir / f"{barcode}")
            if args.fastq:
                path += args.fastq_suffix
                if args.gzip:
                    path += ".gz"
                    handles[barcode] = gzip.open(path, "wt")
                else:
                    handles[barcode] = open(path, "w")
            else:
                path += ".bam"
                handles[barcode] = pysam.AlignmentFile(
                    path, mode="wb", template=bam_file
                )

            fp = handles[barcode]

        if args.fastq:
            read = to_fastq(read)

        fp.write(read)
        barcode_counts[barcode] += 1

    logger.warning(
        f"Number of duplex reads with mismatched barcodes: {num_mismatched_barcodes} - use --verbose to see details"
    )
    logger.warning(
        f"Number of unknown read IDs: {num_orphan_reads} - use --verbose to see details"
    )

    logger.info(f"Total number of reads assessed: {sum(barcode_counts.values())}")

    bam_file.close()

    logger.info("Barcode counts:")
    # sort by barcode
    for barcode, count in sorted(barcode_counts.items(), key=lambda x: x[0]):
        logger.info(f"{barcode}\t{count}")

    # close all the files
    for handle in handles.values():
        handle.close()

    logger.info("Done")


if __name__ == "__main__":
    main()
