import argparse
import pysam


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
        "--output",
        help="Output file path [default: stdout]. Stdout will be a BAM file.",
        default="-",
    )

    args = parser.parse_args()

    # Your code logic goes here
    summary_file = args.summary_file

    mode = "rb"
    if args.bam_file != "-" and args.bam_file.endswith(".sam"):
        mode = "r"
    bam_file = pysam.AlignmentFile(args.bam_file, mode=mode)
    keep_duplex = args.duplex

    mode = "wb"
    if args.output != "-" and args.output.endswith(".sam"):
        mode = "w"
    output_file = pysam.AlignmentFile(args.output, mode=mode, template=bam_file)


if __name__ == "__main__":
    main()
