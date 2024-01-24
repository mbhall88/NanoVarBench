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
        "--outdir",
        required=True,
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
    fixchr_aln = tmpdir / "fixchr.paf"
    cmd = (
        f"minimap2 -x asm5 -t {args.threads} --eqx -c --cs {args.ref} {args.donor} | "
        f"sort -k6,6 -k8,8n > {fixchr_aln}"
    )
    logger.debug("Running minimap2 command: {}", cmd)
    proc = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
    if proc.returncode != 0:
        logger.error("minimap2 failed with error:\n{}", proc.stderr.decode())
        sys.exit(1)

    logger.info("Fixing chromosome orientation with fixchr...")
    cmd = f"fixchr -r {args.ref} -q {args.donor} --dir {tmpdir} -c {fixchr_aln} -F P --contig_size 1000"
    logger.debug("Running fixchr command: {}", cmd)
    proc = subprocess.run(cmd, shell=True, capture_output=True)
    if proc.returncode != 0:
        logger.error("fixchr failed with error:\n{}", proc.stderr.decode())
        sys.exit(1)

    ref = tmpdir / "fixchr.ref.filtered.fa"
    donor = tmpdir / "fixchr.qry.filtered.fa"

    logger.info("Running nucmer...")
    delta = tmpdir / "nucmer.delta"
    delta_filt = tmpdir / "nucmer.1delta"
    coords = tmpdir / "nucmer.coords"
    cmd = f"nucmer -t {args.threads} --maxmatch --prefix={tmpdir}/nucmer --delta={delta} {ref} {donor}"
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
        f"minimap2 -x asm5 -t {args.threads} --eqx -c --cs {ref} {donor} | "
        f"sort -k6,6 -k8,8n > {mm2_aln}"
    )
    logger.debug("Running minimap2 command: {}", cmd)
    proc = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
    if proc.returncode != 0:
        logger.error("minimap2 failed with error:\n{}", proc.stderr.decode())
        sys.exit(1)

    logger.info("Running syri on nucmer alignment...")
    cmd = f"syri --nc {args.threads} --dir {tmpdir} --prefix nucmer. -c {coords} -d {delta_filt} -r {ref} -q {donor}"
    logger.debug("Running syri command: {}", cmd)
    proc = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
    if proc.returncode != 0:
        logger.error("syri failed with error:\n{}", proc.stderr.decode())
        sys.exit(1)
    nucmer_syri = tmpdir / "nucmer.syri.out"
    if not nucmer_syri.exists():
        logger.error(f"syri failed to produce expected output {nucmer_syri}")
        sys.exit(1)

    logger.info("Running syri on minimap2 alignment...")
    cmd = f"syri --nc {args.threads} --dir {tmpdir} --prefix minimap2. -c {mm2_aln} -r {ref} -q {donor} -F P"
    logger.debug("Running syri command: {}", cmd)
    proc = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
    if proc.returncode != 0:
        logger.error("syri failed with error:\n{}", proc.stderr.decode())
        sys.exit(1)
    mm2_syri = tmpdir / "minimap2.syri.out"
    if not mm2_syri.exists():
        logger.error(f"syri failed to produce expected output {mm2_syri}")
        sys.exit(1)

    logger.info("Running plotsr on nucmer alignment...")
    genomes = tmpdir / "genomes.txt"
    with open(genomes, "w") as f:
        print(f"{ref}\treference\tlw:1.5", file=f)
        print(f"{donor}\tdonor\tlw:1.5", file=f)

    nucmer_plot = outdir / "nucmer.plotsr.png"
    cmd = f"plotsr --sr {nucmer_syri} -o {nucmer_plot} --genomes {genomes}"
    logger.debug("Running plotsr command: {}", cmd)
    proc = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
    if proc.returncode != 0:
        logger.error("plotsr failed with error:\n{}", proc.stderr.decode())
        sys.exit(1)

    logger.info("Running plotsr on minimap2 alignment...")
    mm2_plot = outdir / "minimap2.plotsr.png"
    cmd = f"plotsr --sr {mm2_syri} -o {mm2_plot} --genomes {genomes}"
    logger.debug("Running plotsr command: {}", cmd)
    proc = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
    if proc.returncode != 0:
        logger.error("plotsr failed with error:\n{}", proc.stderr.decode())
        sys.exit(1)

    logger.success("Done!")


if __name__ == "__main__":
    main()
