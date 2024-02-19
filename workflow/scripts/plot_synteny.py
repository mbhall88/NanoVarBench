"""This script takes the reference genome and the donor genome, aligns them with nucmer 
and also minimap2, finds synteny and structural rearrangements, for each alignment 
with syri, and then plots each alignment with plotsr.
"""

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

from loguru import logger

LOG_FMT = (
    "[<green>{time:YYYY-MM-DD HH:mm:ss}</green>] <level>{level: <8}</level> | "
    "<level>{message}</level>"
)
PLOTSR_CONFIG = """## COLOURS and transparency for alignments (syntenic, inverted, translocated, and duplicated)
syncol:#CCCCCC
invcol:#FFA500
tracol:#9ACD32
dupcol:#00BBFF
alpha:0.8

## Margins and dimensions:
chrmar:0.1              ## Adjusts the gap between chromosomes and tracks. Higher values leads to more gap
exmar:0.1               ## Extra margin at the top and bottom of plot area

## LEGEND
legend:T                ## To plot legend use T, use F to not plot legend
genlegcol:1            ## Number of columns for genome legend, set -1 for automatic setup
bbox:0,1.01,0.5,0.3     ## [Left edge, bottom edge, width, height]
bbox_v:0,1.1,0.5,0.3    ## For vertical chromosomes (using -v option)
bboxmar:0.5             ## Margin between genome and annotation legends"""


def setup_logging(verbose: bool, quiet: bool) -> None:
    log_lvl = "INFO"
    if verbose:
        log_lvl = "DEBUG"
    elif quiet:
        log_lvl = "ERROR"
    logger.remove()
    logger.add(sys.stderr, level=log_lvl, format=LOG_FMT)


def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--ref",
        required=True,
        help="Path to the reference genome",
    )
    parser.add_argument(
        "--donor",
        required=True,
        help="Path to the donor genome",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        default="./",
        help="Path to the output directory",
    )
    parser.add_argument(
        "-t",
        "--threads",
        default=1,
        help="Number of threads to use",
    )
    parser.add_argument(
        "--tmpdir",
        help="Path to the temporary directory",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Print debug information",
    )
    parser.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Only print errors",
    )
    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    setup_logging(args.verbose, args.quiet)

    # Create output directory
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    # Create temporary directory
    if args.tmpdir:
        tmpdir = Path(args.tmpdir).resolve()
        tmpdir.mkdir(parents=True, exist_ok=True)
    else:
        tmpdir = Path(tempfile.mkdtemp())

    # align with minimap2 and then fix chromosome orientation with fixchr
    logger.info("Fix chromosome orientation...")
    chroder_aln = tmpdir / "chroder.sam"
    cmd = (
        f"minimap2 --eqx --cs -ax asm5 -t {args.threads} {args.ref} {args.donor} > {chroder_aln}"
        # f"minimap2 -x asm5 -t {args.threads} --eqx -c --cs {args.ref} {args.donor} | "
        # f"sort -k6,6 -k8,8n > {fixchr_aln}"
    )
    logger.debug("Running minimap2 command: {}", cmd)
    proc = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
    if proc.returncode != 0:
        logger.error("minimap2 failed with error:\n{}", proc.stderr.decode())
        sys.exit(1)

    logger.info("Generating psuedo-chromosome(s) for donor...")
    prefix = tmpdir / "chroder"
    cmd = f"chroder -F S -o {prefix} {chroder_aln} {args.ref} {args.donor}"
    logger.debug("Running chroder command: {}", cmd)
    proc = subprocess.run(
        cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
    )
    if proc.returncode != 0:
        logger.error("chroder failed with error:\n{}", proc.stdout.decode())
        sys.exit(1)

    # ref = tmpdir / f"{prefix}fixchr.ref.filtered.fa"
    ref = args.ref
    donor = prefix.with_suffix(".qry.fasta")

    # filter ref and donor so they have the same chromosomes
    logger.debug("Filtering reference and donor genomes...")
    ref_filtered = tmpdir / "ref.filtered.fa"
    donor_filtered = tmpdir / "donor.filtered.fa"
    ref_chrs = set()
    donor_chrs = set()
    with open(ref) as f:
        for line in f:
            if line.startswith(">"):
                ref_chrs.add(line.split()[0][1:])
    with open(donor) as f:
        for line in f:
            if line.startswith(">"):
                donor_chrs.add(line.split()[0][1:])

    common_chrs = ref_chrs & donor_chrs
    if not common_chrs:
        logger.error("No common chromosomes found between reference and donor")
        sys.exit(1)
    else:
        logger.debug("Common chromosomes: {}", common_chrs)

    with open(ref_filtered, "w") as out:
        with open(ref) as f:
            keep_seq = False
            for line in f:
                if line.startswith(">"):
                    if line.split()[0][1:] in common_chrs:
                        keep_seq = True
                        out.write(line)
                    else:
                        keep_seq = False
                elif keep_seq:
                    out.write(line)

    with open(donor_filtered, "w") as out:
        with open(donor) as f:
            keep_seq = False
            for line in f:
                if line.startswith(">"):
                    if line.split()[0][1:] in common_chrs:
                        keep_seq = True
                        out.write(line)
                    else:
                        keep_seq = False
                elif keep_seq:
                    out.write(line)

    logger.info("Running nucmer...")
    delta = tmpdir / "nucmer.delta"
    delta_filt = tmpdir / "nucmer.1delta"
    coords = tmpdir / "nucmer.coords"
    cmd = f"nucmer -t {args.threads} --maxmatch --prefix={tmpdir}/nucmer {ref_filtered} {donor_filtered}"
    logger.debug("Running nucmer command: {}", cmd)
    proc = subprocess.run(cmd, shell=True, capture_output=True)
    if proc.returncode != 0:
        logger.error("nucmer failed with error:\n{}", proc.stderr.decode())
        sys.exit(1)

    logger.info("Running delta-filter...")
    cmd = f"delta-filter -1 {delta} > {delta_filt}"
    logger.debug("Running delta-filter command: {}", cmd)
    proc = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
    if proc.returncode != 0:
        logger.error("delta-filter failed with error:\n{}", proc.stderr.decode())
        sys.exit(1)

    logger.info("Running show-coords...")
    cmd = f"show-coords -THrd {delta_filt} > {coords}"
    logger.debug("Running show-coords command: {}", cmd)
    proc = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
    if proc.returncode != 0:
        logger.error("show-coords failed with error:\n{}", proc.stderr.decode())
        sys.exit(1)

    logger.info("Running minimap2...")
    mm2_aln = tmpdir / "minimap2.paf"
    cmd = (
        f"minimap2 -x asm5 -t {args.threads} --eqx -c --cs {ref_filtered} {donor_filtered} | "
        f"sort -k6,6 -k8,8n > {mm2_aln}"
    )
    logger.debug("Running minimap2 command: {}", cmd)
    proc = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
    if proc.returncode != 0:
        logger.error("minimap2 failed with error:\n{}", proc.stderr.decode())
        sys.exit(1)

    logger.info("Running syri on nucmer alignment...")
    cmd = f"syri --nc {args.threads} --dir {tmpdir} --prefix nucmer. -c {coords} -d {delta_filt} -r {ref_filtered} -q {donor_filtered}"
    logger.debug("Running syri command: {}", cmd)
    proc = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
    if proc.returncode != 0:
        logger.error("syri failed with error:\n{}", proc.stderr.decode())
        sys.exit(1)
    else:
        syri_summary = tmpdir / "nucmer.syri.summary"
        logger.debug("nucmer syri summary:\n{}", syri_summary.read_text())

    nucmer_syri = tmpdir / "nucmer.syri.out"
    if not nucmer_syri.exists():
        logger.error(f"syri failed to produce expected output {nucmer_syri}")
        sys.exit(1)

    logger.info("Running syri on minimap2 alignment...")
    cmd = f"syri --nc {args.threads} --dir {tmpdir} --prefix minimap2. -c {mm2_aln} -r {ref_filtered} -q {donor_filtered} -F P"
    logger.debug("Running syri command: {}", cmd)
    proc = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)

    if proc.returncode != 0:
        logger.error("syri failed with error:\n{}", proc.stderr.decode())
        sys.exit(1)
    else:
        syri_summary = tmpdir / "minimap2.syri.summary"
        logger.debug("minimap2 syri summary:\n{}", syri_summary.read_text())

    mm2_syri = tmpdir / "minimap2.syri.out"
    if not mm2_syri.exists():
        logger.error(f"syri failed to produce expected output {mm2_syri}")
        sys.exit(1)

    logger.info("Running plotsr on nucmer alignment...")
    config = tmpdir / "plotsr.config"
    config.write_text(PLOTSR_CONFIG)
    plotsr_opts = f"--cfg {config}"
    genomes = tmpdir / "genomes.txt"
    with open(genomes, "w") as f:
        print(f"{ref_filtered}\treference\tlw:1.5", file=f)
        print(f"{donor_filtered}\tdonor\tlw:1.5", file=f)

    nucmer_plot = outdir / "nucmer.plotsr.png"
    cmd = (
        f"plotsr {plotsr_opts} --sr {nucmer_syri} -o {nucmer_plot} --genomes {genomes}"
    )
    logger.debug("Running plotsr command: {}", cmd)
    proc = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
    if proc.returncode != 0:
        logger.error("plotsr failed with error:\n{}", proc.stderr.decode())
        sys.exit(1)
    else:
        logger.success("Nucmer plotsr figure saved to {}", nucmer_plot)

    logger.info("Running plotsr on minimap2 alignment...")
    mm2_plot = outdir / "minimap2.plotsr.png"
    cmd = f"plotsr {plotsr_opts} --sr {mm2_syri} -o {mm2_plot} --genomes {genomes}"
    logger.debug("Running plotsr command: {}", cmd)
    proc = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
    if proc.returncode != 0:
        logger.error("plotsr failed with error:\n{}", proc.stderr.decode())
        sys.exit(1)
    else:
        logger.success("Minimap2 plotsr figure saved to {}", mm2_plot)

    logger.success("Done!")


if __name__ == "__main__":
    main()
