"""This script creates a truth VCF for a reference genome. Variants are taken from a 
donor genome of the given species, with a given ANI from the reference 
genome. The variants are generated using minimap2 and dnadiff, and then merged, 
normalised, and filtered using bcftools. The final VCF is then used to create a mutant 
reference genome using bcftools consensus.
"""

import argparse
import atexit
import gzip
import json
import operator
import random
import re
import shutil
import subprocess
import sys
import tempfile
from collections import OrderedDict, defaultdict
from pathlib import Path
from typing import Dict, List, Tuple
from uuid import uuid4

import pyfastaq
from intervaltree import Interval, IntervalTree
from loguru import logger

MIN_BCFTOOLS_VERSION = "1.18"
LOG_FMT = (
    "[<green>{time:YYYY-MM-DD HH:mm:ss}</green>] <level>{level: <8}</level> | "
    "<level>{message}</level>"
)
FASTA_REGEX = re.compile(r"^.*(\.fa|\.fasta|\.fna)(\.gz)?$")
ACC_REGEX = re.compile(r"(?P<acc>GC[AF]_\d+\.\d+)")
SNP = 1
DEL = 2
INS = 3

var_types = {
    1: "SNP",
    2: "DEL",
    3: "INS",
}


def setup_logging(verbose: int = 0, quiet: bool = False) -> None:
    log_lvl = "INFO"
    if verbose == 1:
        log_lvl = "DEBUG"
    elif verbose > 1:
        log_lvl = "TRACE"
    elif quiet:
        log_lvl = "ERROR"
    logger.remove()
    logger.add(sys.stderr, level=log_lvl, format=LOG_FMT)


def parse_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=__doc__,
    )
    parser.add_argument("reference_genome", type=str, help="Reference genome")
    parser.add_argument(
        "-s",
        "--species",
        type=str,
        required=True,
        help='Species - e.g., "s__Escherichia coli" (GTDB) or "Escherichia coli" (NCBI)',
    )
    parser.add_argument(
        "-a",
        "--ani",
        type=float,
        default=99.50,
        help="Pick genome with ANI closest to this value",
    )
    parser.add_argument(
        "-m",
        "--min-ani",
        type=float,
        default=0.0,
        help="Minimum ANI to consider",
    )
    parser.add_argument(
        "-M",
        "--max-ani",
        type=float,
        default=100.0,
        help="Maximum ANI to consider",
    )
    parser.add_argument(
        "-C",
        "--min-completeness",
        type=float,
        default=0.0,
        help="Minimum assembly (CheckM) completeness (%%) to consider",
    )
    parser.add_argument(
        "-P",
        "--min-completeness-percentile",
        type=float,
        default=0.0,
        help="Minimum assembly (CheckM) completeness percentile (%%) to consider",
    )
    parser.add_argument(
        "-X",
        "--max-contamination",
        type=float,
        default=100.0,
        help="Maximum assembly (CheckM) contamination (%%) to consider",
    )
    parser.add_argument(
        "-c",
        "--no-cleanup",
        action="store_true",
        help="Don't delete temporary files",
    )
    parser.add_argument(
        "-t", "--threads", type=int, default=1, help="Number of threads to use"
    )
    parser.add_argument(
        "-l",
        "--asmlvl",
        default="complete genome",
        help="Assembly level. (comma-separated entries, empty ('') for all) [complete genome, chromosome, scaffold, contig]",
    )
    parser.add_argument(
        "-T",
        "--taxonomy",
        default="gtdb",
        choices=["gtdb", "ncbi"],
        help="Taxonomy database to use",
    )
    parser.add_argument(
        "-D",
        "--database",
        default="refseq",
        help="Database to use (comma-separated entries) [refseq, genbank]",
    )
    parser.add_argument(
        "-A",
        "--max-asm",
        type=int,
        default=10_000,
        help="Maximum number of assemblies to download. 0 for all",
    )
    parser.add_argument(
        "-o", "--outdir", type=str, default="./mutref", help="Output directory"
    )
    parser.add_argument("--tmpdir", type=str, help="Temporary directory")
    parser.add_argument(
        "-f",
        "--existing",
        help="Directory containing existing genomes (skips downloading genomes)",
    )
    parser.add_argument(
        "-I", "--max-indel", type=int, default=50, help="Maximum indel size"
    )
    parser.add_argument(
        "-O",
        "--remove-overlaps",
        action="store_true",
        help="Remove valid variants whose ALT start and end overlap (see https://github.com/samtools/bcftools/issues/2082)",
    )
    parser.add_argument(
        "-F",
        "--force",
        action="store_true",
        help="Force overwrite existing output directory",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="Increase verbosity level. Can be used multiple times",
    )
    parser.add_argument("-q", "--quiet", action="store_true", help="Quiet logging")
    return parser.parse_args()


# add function to test if command is executable
def is_executable(cmd: str) -> bool:
    path = shutil.which(cmd)
    return path is not None


def check_bcftools_version(min_version: str) -> None:
    bcftools_version = subprocess.check_output(["bcftools", "--version"]).decode()
    bcftools_version = re.search(r"bcftools (\d+\.\d+)", bcftools_version).group(1)

    if bcftools_version is None:
        logger.error("Could not determine bcftools version")
        sys.exit(1)

    if bcftools_version < min_version:
        logger.error(
            f"bcftools version {min_version} or greater is required. You have {bcftools_version}"
        )
        sys.exit(1)


def download_genomes(
    tmpdir: str,
    database: str,
    taxonomy: str,
    species: str,
    threads: int,
    asmlvl: str,
    max_asm: int = 0,
):
    args = [
        "genome_updater.sh",
        "-d",
        database,
        "-g",
        "bacteria",
        "-f",
        "genomic.fna.gz",
        "-o",
        tmpdir,
        "-M",
        taxonomy,
        "-T",
        species,
        "-A",
        f"species:{max_asm}",
        "-m",
        "-a",
        "-t",
        str(threads),
        "-l",
        asmlvl,
    ]
    logger.debug(f"Running genome_updater.sh with args: {' '.join(args)}")
    proc = subprocess.run(
        args,
        stderr=subprocess.STDOUT,
        stdout=subprocess.PIPE,
        text=True,
    )

    if proc.returncode != 0:
        logger.error("Error downloading genomes")
        logger.error(proc.stdout)
        sys.exit(1)
    else:
        logger.trace(f"genome_updater.sh output:\n{proc.stdout}")


def sketch_genomes(genomes_fofn: str, sketch_dir: str, threads: int) -> None:
    args = ["skani", "sketch", "-o", sketch_dir, "-t", str(threads), "-l", genomes_fofn]
    logger.debug(f"Running skani sketch with args: {' '.join(args)}")
    proc = subprocess.run(
        args,
        stderr=subprocess.STDOUT,
        stdout=subprocess.PIPE,
        text=True,
    )

    if proc.returncode != 0:
        logger.error("Error sketching genomes")
        logger.error(proc.stdout)
        sys.exit(1)
    else:
        logger.trace(f"skani sketch output:\n{proc.stdout}")


def calculate_ANI(
    sketch_dir: str, query_genome: str, threads: int, output: str
) -> None:
    args = [
        "skani",
        "search",
        "-o",
        output,
        "-t",
        str(threads),
        "-d",
        sketch_dir,
        "-q",
        query_genome,
    ]
    logger.debug(f"Running skani with args: {' '.join(args)}")
    proc = subprocess.run(
        args,
        stderr=subprocess.STDOUT,
        stdout=subprocess.PIPE,
        text=True,
    )

    if proc.returncode != 0:
        logger.error("Error running skani")
        logger.error(proc.stdout)
        sys.exit(1)
    else:
        logger.trace(f"skani output:\n{proc.stdout}")


def argclosest(array: List[Tuple[float, int]], value: float) -> Tuple[float, int]:
    """Takes a list of tuples (distance, index) and returns the index of the distance closest to the given value"""
    argmin = min(range(len(array)), key=lambda i: abs(array[i][0] - value))
    return array[argmin]


def is_file_gzipped(file: str) -> bool:
    """Checks if a file is gzipped"""
    return open(file, "rb").read(2) == b"\x1f\x8b"


def generate_minimap2_variants(
    reference_genome: str, donor_genome: str, vcf: str, threads: int
) -> None:
    cmd = (
        f"minimap2 -x asm5 -c --cs {reference_genome} {donor_genome} |"
        "sort -k6,6 -k8,8n |"
        f"paftools.js call -f {reference_genome} -l50 -L50 - |"
        f"bcftools sort -o {vcf} --write-index -"
    )
    logger.debug(f"Running minimap2 with args: {cmd}")

    proc = subprocess.run(
        cmd,
        shell=True,
        capture_output=True,
        text=True,
    )

    if proc.returncode != 0:
        logger.error("Error running minimap2")
        logger.error(proc.stderr)
        sys.exit(1)
    else:
        logger.trace(f"minimap2 output:\n{proc.stderr}")


def generate_dnadiff_variants(
    tmpdir: Path, reference_genome: str, donor_genome: str, vcf: str, threads: int
) -> None:
    snps_file = tmpdir / "dnadiff.snps"
    run_dnadiff(
        tmpdir,
        reference_genome,
        donor_genome,
        snps_file,
        threads=threads,
    )
    raw_vcf = tmpdir / "dnadiff.vcf"
    snps_file_to_vcf(snps_file, reference_genome, raw_vcf)

    # compress, sort and index raw VCF
    cmd = f"bcftools sort -o {vcf} --write-index {raw_vcf}"
    logger.debug(f"Running bcftools sort with args: {cmd}")
    proc = subprocess.run(
        cmd,
        shell=True,
        capture_output=True,
        text=True,
    )

    if proc.returncode != 0:
        logger.error("Error running bcftools sort")
        logger.error(proc.stderr)
        sys.exit(1)
    else:
        logger.trace(f"bcftools sort output:\n{proc.stderr}")


def run_dnadiff(
    tmpdir: Path,
    reference_genome: str,
    donor_genome: str,
    snps_file: str,
    threads=1,
):
    delta = tmpdir / "nucmer.delta"
    delta_1 = tmpdir / "nucmer.1delta"
    if delta.exists():
        delta.unlink()
    if delta_1.exists():
        delta_1.unlink()

    nucmer_cmd = f"nucmer --maxmatch --threads {threads} --delta {delta} {donor_genome} {reference_genome}"

    logger.debug(f"Running nucmer with args: {nucmer_cmd}")
    nucmer_proc = subprocess.run(nucmer_cmd, shell=True, capture_output=True, text=True)

    if nucmer_proc.returncode != 0:
        logger.error("Error running nucmer")
        logger.error(nucmer_proc.stderr)
        sys.exit(1)
    else:
        logger.trace(f"nucmer output:\n{nucmer_proc.stderr}")

    filter_cmd = f"delta-filter -1 {delta} > {delta_1}"

    logger.debug(f"Running delta-filter with args: {filter_cmd}")
    filter_proc = subprocess.run(filter_cmd, shell=True, capture_output=True, text=True)

    if filter_proc.returncode != 0:
        logger.error("Error running delta-filter")
        logger.error(filter_proc.stderr)
        sys.exit(1)
    else:
        logger.trace(f"delta-filter output:\n{filter_proc.stderr}")

    snps_cmd = f"show-snps -rlTHC {delta_1} > {snps_file}"

    logger.debug(f"Running show-snps with args: {snps_cmd}")
    snps_proc = subprocess.run(snps_cmd, shell=True, capture_output=True, text=True)

    if snps_proc.returncode != 0:
        logger.error("Error running show-snps")
        logger.error(snps_proc.stderr)
        sys.exit(1)
    else:
        logger.trace(f"show-snps output:\n{snps_proc.stderr}")


def snps_file_to_vcf(snps_file: str, reference_genome: str, vcf: str):
    vcf_records = {}
    variants = get_all_mummer_snps(snps_file)
    ref_seqs = file_to_dict_of_seqs(reference_genome)

    for variant in variants:
        # If the variant is reversed, it means that either the ref or query had to be
        # reverse complemented when aligned by mummer. Need to do the appropriate
        # reverse (complement) fixes so the VCF has the correct REF and ALT sequences
        if variant.reverse:
            qry_seq = pyfastaq.sequences.Fasta("x", variant.qry_base)
            qry_seq.revcomp()
            variant.qry_base = "".join(reversed(qry_seq.seq))
            ref_seq = pyfastaq.sequences.Fasta("x", variant.ref_base)
            ref_seq.revcomp()
            variant.ref_base = ref_seq.seq

        if variant.var_type == SNP:
            new_record = VcfRecord(
                "\t".join(
                    [
                        variant.qry_name,
                        str(variant.qry_start + 1),
                        ".",
                        variant.qry_base,
                        variant.ref_base,
                        ".",
                        ".",
                        "VTYPE=DNADIFF_SNP",
                        "GT",
                        "1/1",
                    ]
                )
            )
        elif variant.var_type == DEL:
            # The query has sequence missing, compared to the reference.
            new_record = VcfRecord(
                "\t".join(
                    [
                        variant.qry_name,
                        str(variant.qry_start + 1),
                        ".",
                        ref_seqs[variant.qry_name][variant.qry_start],
                        ref_seqs[variant.qry_name][variant.qry_start]
                        + variant.ref_base,
                        ".",
                        ".",
                        "VTYPE=DNADIFF_INS",
                        "GT",
                        "1/1",
                    ]
                )
            )
        elif variant.var_type == INS:
            # The ref has sequence missing, compared to the query.
            new_record = VcfRecord(
                "\t".join(
                    [
                        variant.qry_name,
                        str(variant.qry_start),
                        ".",
                        ref_seqs[variant.qry_name][variant.qry_start - 1]
                        + variant.qry_base,
                        ref_seqs[variant.qry_name][variant.qry_start - 1],
                        ".",
                        ".",
                        "VTYPE=DNADIFF_DEL",
                        "GT",
                        "1/1",
                    ]
                )
            )
        else:
            raise Exception("Unknown variant type: " + str(variant))

        new_record.INFO["QNAME"] = variant.ref_name
        new_record.INFO["QSTART"] = str(variant.ref_start + 1)
        new_record.INFO["QSTRAND"] = "-" if variant.reverse else "+"

        assert (
            new_record.REF
            == ref_seqs[new_record.CHROM][
                new_record.POS : new_record.POS + len(new_record.REF)
            ]
        )

        if new_record.CHROM not in vcf_records:
            vcf_records[new_record.CHROM] = []

        vcf_records[new_record.CHROM].append(new_record)

    for vcf_list in vcf_records.values():
        vcf_list.sort(key=operator.attrgetter("POS"))

    with open(vcf, "w") as f:
        print("##fileformat=VCFv4.2", file=f)
        print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">', file=f)
        print(
            '##INFO=<ID=VTYPE,Number=1,Type=String,Description="Type of variant">',
            file=f,
        )
        print('##INFO=<ID=QNAME,Number=1,Type=String,Description="Query name">', file=f)
        print(
            '##INFO=<ID=QSTART,Number=1,Type=Integer,Description="Query start">',
            file=f,
        )
        print(
            '##INFO=<ID=QSTRAND,Number=1,Type=String,Description="Query strand">',
            file=f,
        )
        for seq in ref_seqs.values():
            print(f"##contig=<ID={seq.id},length={len(seq)}>", file=f)
        print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample", file=f)

        for key, vcf_list in sorted(vcf_records.items()):
            for record in vcf_list:
                print(record, file=f)


def file_to_dict_of_seqs(infile: str):
    """Given a file of sequences, returns a dictionary of
    sequence name -> pyfastaq.sequences.Fasta.
    Anything after first whitespace is removed from the names"""
    seqs = {}
    pyfastaq.tasks.file_to_dict(infile, seqs)
    seqs = {k.split()[0]: v for k, v in seqs.items()}
    for seq in seqs.values():
        seq.id = seq.id.split()[0]
        seq.seq = seq.seq.upper()
    return seqs


def get_all_mummer_snps(fname: str):
    variants = []
    fr = snps_reader(fname)
    for nucmer_snp in fr:
        if len(variants) == 0 or not variants[-1].update_indel(nucmer_snp):
            variants.append(Variant(nucmer_snp))

    return variants


def snps_reader(fname: str):
    fd = open(fname, "r")
    for line in fd:
        if line.startswith("[") or ("\t" not in line):
            continue

        yield Snp(line)

    fd.close()


class Snp:
    def __init__(self, line: str):
        # Without the -C option to show-snps, looks like this:
        # [P1] [SUB] [SUB]  [P2] [BUFF] [DIST] [LEN R] [LEN Q] [FRM]   [TAGS]
        # 187  A     C      269  187    187    654     853     1       1   ref_name  qry_name

        # With the -C option to show-snps, looks like this:
        # [P1] [SUB] [SUB] [P2] [BUFF] [DIST] [R] [Q] [LEN R] [LEN Q] [FRM] [TAGS]
        # 187  A     C     269  187    187    0   0   654     853     1     1   ref_name  qry_name
        try:
            fields = line.rstrip().split("\t")
            self.ref_pos = int(fields[0]) - 1
            self.ref_base = fields[1]
            self.qry_base = fields[2]
            self.qry_pos = int(fields[3]) - 1
            self.ref_length = int(fields[-6])
            self.qry_length = int(fields[-5])
            self.reverse = {"1": False, "-1": True}[fields[-3]]
            self.ref_name = fields[-2]
            self.qry_name = fields[-1]
        except Exception:
            raise ValueError(
                "Error constructing Snp from mummer show-snps output at this line:\n"
                + line
            )

    def __eq__(self, other):
        return type(other) is type(self) and self.__dict__ == other.__dict__

    def __str__(self):
        return "\t".join(
            [
                str(x)
                for x in [
                    self.ref_pos + 1,
                    self.ref_base,
                    self.qry_base,
                    self.qry_pos + 1,
                    self.ref_length,
                    self.qry_length,
                    "-1" if self.reverse else "1",
                    self.ref_name,
                    self.qry_name,
                ]
            ]
        )


class Variant:
    def __init__(self, snp: Snp):
        """Create a Variant object from a pymummer.snp.Snp object"""
        if snp.ref_base == ".":
            self.var_type = INS
            self.qry_base = snp.qry_base
            self.ref_base = "."
        elif snp.qry_base == ".":
            self.var_type = DEL
            self.qry_base = "."
            self.ref_base = snp.ref_base
        elif "." not in [snp.ref_base, snp.qry_base]:
            self.var_type = SNP
            self.ref_base = snp.ref_base
            self.qry_base = snp.qry_base
        else:
            raise ValueError("Error constructing Variant from Snp: " + str(snp))

        self.ref_start = snp.ref_pos
        self.ref_end = snp.ref_pos
        self.ref_length = snp.ref_length
        self.ref_name = snp.ref_name
        self.qry_start = snp.qry_pos
        self.qry_end = snp.qry_pos
        self.qry_length = snp.qry_length
        self.qry_name = snp.qry_name
        self.reverse = snp.reverse

    def __eq__(self, other):
        return type(other) is type(self) and self.__dict__ == other.__dict__

    def __str__(self):
        return "\t".join(
            [
                str(self.ref_start + 1),
                str(self.ref_end + 1),
                str(self.ref_length),
                str(self.ref_name),
                self.ref_base,
                str(self.qry_start + 1),
                str(self.qry_end + 1),
                str(self.qry_length),
                str(self.qry_name),
                self.qry_base,
                "-1" if self.reverse else "1",
            ]
        )

    def update_indel(self, nucmer_snp: Snp) -> bool:
        """Indels are reported over multiple lines, 1 base insertion or deletion per line. This method extends the current variant by 1 base if it's an indel and adjacent to the new SNP and returns True. If the current variant is a SNP, does nothing and returns False"""
        new_variant = Variant(nucmer_snp)
        if (
            self.var_type not in [INS, DEL]
            or self.var_type != new_variant.var_type
            or self.qry_name != new_variant.qry_name
            or self.ref_name != new_variant.ref_name
            or self.reverse != new_variant.reverse
        ):
            return False
        if (
            self.var_type == INS
            and self.ref_start == new_variant.ref_start
            and self.qry_end + 1 == new_variant.qry_start
        ):
            self.qry_base += new_variant.qry_base
            self.qry_end += 1
            return True
        if (
            self.var_type == DEL
            and self.qry_start == new_variant.qry_start
            and self.ref_end + 1 == new_variant.ref_start
        ):
            self.ref_base += new_variant.ref_base
            self.ref_end += 1
            return True

        return False


class VcfRecord:
    def __init__(self, line):
        """Constructs VcfRecord from a line of a VCF file.
        Assumes only one sample in the file"""
        assert not line.startswith("#")
        fields = line.rstrip().split("\t")
        # #CHROM          POS     ID                      REF     ALT     QUAL    FILTER  INFO                            FORMAT          2.2.2.1
        # NC_000962.3     1977    UNION_BC_k31_var_120    A       G       .       PASS    KMER=31;SVLEN=0;SVTYPE=SNP      GT:COV:GT_CONF  1/1:0,52:39.80
        try:
            self.CHROM = fields[0]
            self.POS = int(fields[1]) - 1
            self.ID = fields[2]
            self.REF = fields[3]
            self.ALT = fields[4].split(",")
            self.QUAL = fields[5]
            self.FILTER = set() if fields[6] == "." else set(fields[6].split(";"))
            INFO = fields[7]
        except Exception:
            raise Exception("Error reading line of vcf file:" + line)

        if self.POS < 0:
            raise ValueError(
                f"POS value {self.POS + 1}, which is less than 1. Cannot continue. Line of VCF file:\n{line}"
            )

        try:
            self.QUAL = float(self.QUAL)
        except Exception:
            self.QUAL = None

        self.INFO = {}
        if INFO != ".":
            info_fields = INFO.split(";")
            for field in info_fields:
                if "=" in field:
                    key, value = field.split("=")
                    self.INFO[key] = value
                else:
                    self.INFO[field] = None

        if len(fields) == 10:
            format_keys = fields[8].split(":")
            format_vals = fields[9].split(":")
            self.FORMAT = OrderedDict(zip(format_keys, format_vals))
            if "GT" in self.FORMAT:
                self.FORMAT.move_to_end("GT", last=False)
        else:
            self.FORMAT = OrderedDict()

    def __eq__(self, other):
        return type(other) is type(self) and self.__dict__ == other.__dict__

    def __repr__(self):
        if len(self.INFO) == 0:
            info_string = "."
        else:
            info_fields = []
            for x in sorted(self.INFO):
                if self.INFO[x] is None:
                    info_fields.append(x)
                else:
                    info_fields.append(x + "=" + self.INFO[x])
            info_string = ";".join(info_fields)

        fields = [
            self.CHROM,
            str(self.POS + 1),
            self.ID,
            self.REF,
            ",".join(self.ALT),
            "." if self.QUAL is None else str(self.QUAL),
            "." if len(self.FILTER) == 0 else ";".join(sorted(list(self.FILTER))),
            info_string,
        ]

        if len(self.FORMAT) > 0:
            fields.append(":".join(self.FORMAT.keys()))
            fields.append(":".join(self.FORMAT.values()))

        return "\t".join(fields)


def merge_and_filter_variants(
    tmpdir: Path,
    vcfs: List[str],
    reference_genome: str,
    output: str,
    max_indel: int,
    remove_valid_overlaps: bool = False,
) -> None:
    tmpvcf = tmpdir / "merged.vcf"
    cmd = (
        f"bcftools concat -Da {' '.join(vcfs)} | "
        f"bcftools norm -f {reference_genome} -a -c e -m - |"
        "bcftools norm -aD |"
        "bcftools +remove-overlaps - |"
        "bcftools annotate --remove QUAL,INFO/VTYPE,INFO/QNAME,INFO/QSTRAND,INFO/QSTART - |"
        "bcftools +setGT - -- -t a -n c:1 |"  # make genotypes haploid ALT
        f"bcftools filter -e 'abs(ILEN)>{max_indel}' -o {tmpvcf}"
    )
    logger.debug(f"Running bcftools pipeline: {cmd}")
    proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    if proc.returncode != 0:
        logger.error(f"bcftools pipeline failed with return code {proc.returncode}")
        logger.error(proc.stderr)
        sys.exit(1)
    else:
        logger.trace(f"bcftools pipeline output:\n{proc.stderr}")

    logger.debug("Adding UIDs to VCF records")

    tmp_uid_vcf = tmpdir / "uid.vcf"
    # add UIDs to every VCF record
    with open(tmp_uid_vcf, "w") as vcf_out, open(tmpvcf, "r") as vcf_in:
        new_ids = set()
        for i, line in enumerate(vcf_in):
            if line.startswith("#"):
                vcf_out.write(line)
            else:
                fields = line.split("\t")
                current_id = fields[2]
                if current_id != ".":
                    vcf_out.write(line)
                else:
                    new_id = shortuuid(seed=i)
                    while new_id in new_ids:
                        new_id = shortuuid(seed=i)
                    new_ids.add(new_id)
                    fields[2] = new_id
                    vcf_out.write("\t".join(fields))

    if remove_valid_overlaps:
        logger.debug("Removing (valid) overlapping variants...")
        rm_overlap_vcf = tmpdir / "rm_overlap.vcf"
        remove_overlaps(tmp_uid_vcf, rm_overlap_vcf)
    else:
        rm_overlap_vcf = tmp_uid_vcf

    cmd = f"bcftools view --write-index -o {output} {rm_overlap_vcf}"
    logger.debug(f"Running bcftools view: {cmd}")
    proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    if proc.returncode != 0:
        logger.error(f"bcftools view failed with return code {proc.returncode}")
        logger.error(proc.stderr)
        sys.exit(1)
    else:
        logger.trace(f"bcftools view output:\n{proc.stderr}")


def remove_overlaps(input: str, output: str) -> None:
    """Removes overlapping variants from a VCF file"""
    intervals = defaultdict(IntervalTree)
    n_records = 0
    with open(input, "r") as vcf_in:
        for line in vcf_in:
            if line.startswith("#"):
                continue
            n_records += 1
            fields = line.split("\t")
            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]
            vid = fields[2]

            if "," in alt:
                logger.error(
                    f"Cannot remove overlapping variants from VCF with multiple ALTs at line:\n{line}"
                )
                sys.exit(1)
            if len(ref) == len(alt):
                # SNP
                iv = Interval(pos, pos + len(ref), vid)
            elif len(ref) > len(alt):
                # deletion
                iv = Interval(pos, pos + len(ref), vid)
            else:
                # insertion
                iv = Interval(pos, pos + len(alt), vid)

            intervals[chrom].add(iv)

    ids_to_keep = non_overlapping_intervals(intervals)
    with open(output, "w") as vcf_out, open(input, "r") as vcf_in:
        for line in vcf_in:
            if line.startswith("#"):
                vcf_out.write(line)
            else:
                vid = line.split("\t")[2]
                if vid in ids_to_keep:
                    vcf_out.write(line)

    logger.debug(f"Removed {n_records - len(ids_to_keep)} overlapping variants")


def non_overlapping_intervals(tree: Dict[str, IntervalTree]) -> IntervalTree:
    keep = set()
    for tree in tree.values():
        for iv in tree:
            vid = iv.data
            overlaps = list(tree.overlap(iv))
            if len(overlaps) == 1:
                assert overlaps[0].data == vid
                keep.add(vid)
            elif len(overlaps) == 0:
                logger.error(
                    f"No overlapping variants found for {vid} - not even with itself"
                )
                sys.exit(1)

    return keep


def shortuuid(length: int = 8, seed=None) -> str:
    random.seed(seed)
    s = str(uuid4()) + str(uuid4())
    s = s.replace("-", "")
    return "".join(random.sample(s, length))


def vcfstats(vcf: str, output: str, log_stats: bool = True) -> None:
    cmd = f"paftools.js vcfstat {vcf}"
    if log_stats:
        cmd += f" | tee {output}"
        proc = subprocess.run(
            cmd,
            shell=True,
            text=True,
            stderr=subprocess.STDOUT,
            stdout=subprocess.PIPE,
        )
        if proc.returncode != 0:
            logger.error(
                f"paftools.js vcfstat failed with return code {proc.returncode}"
            )
            logger.error(proc.stdout)
            sys.exit(1)
        else:
            logger.info("VCF stats:")
            for line in proc.stdout.split("\n"):
                if not line:
                    continue
                logger.info(line)
    else:
        cmd += f" > {output}"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if proc.returncode != 0:
            logger.error(
                f"paftools.js vcfstat failed with return code {proc.returncode}"
            )
            logger.error(proc.stderr)
            sys.exit(1)


def create_mutant_reference(reference_genome: str, variants: str, output: str) -> None:
    cmd = f"bcftools consensus -f {reference_genome} {variants} > {output}"
    logger.debug(f"Running bcftools consensus: {cmd}")
    proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    if proc.returncode != 0:
        logger.error(f"bcftools consensus failed with return code {proc.returncode}")
        logger.error(proc.stderr)
        sys.exit(1)
    else:
        log = proc.stderr
        if "overlaps with another variant" in log:
            logger.error(
                "bcftools consensus warned there are overlapping variants. This should not happen."
            )
            logger.error(log)
            sys.exit(1)
        logger.trace(f"bcftools consensus output:\n{proc.stderr}")


def inverse_vcf(invcf: str, outvcf: str) -> None:
    """This function takes a VCF that was used to generate a consensus sequence with bcftools
    consensus and provides the inverse - i.e., the VCF that would turn the consensus
    sequence back into the original reference sequence.
    """
    with open(outvcf, "w") as f_out, gzip.open(invcf, "rt") as f_in:
        offset_counter = defaultdict(int)
        for line in f_in:
            if line.startswith("#"):
                print(line.strip(), file=f_out)
                continue

            fields = line.strip().split("\t")
            ref = fields[3]
            alt = fields[4]
            pos = int(fields[1])
            chrom = fields[0]
            offset = offset_counter[chrom]
            new_alt = ref
            new_ref = alt
            new_pos = pos + offset
            record_offset = len(alt) - len(ref)
            offset_counter[chrom] += record_offset
            fields[1] = str(new_pos)
            fields[3] = new_ref
            fields[4] = new_alt
            print("\t".join(fields), file=f_out)


def fasta_to_dict(fasta: str) -> Dict[str, str]:
    """Read a fasta file and return a dictionary of sequences"""
    seqs = {}
    with open(fasta, "r") as f:
        for line in f:
            if line.startswith(">"):
                name = line.strip().split()[0][1:]
                seqs[name] = ""
            else:
                seqs[name] += line.strip()

    return seqs


def fastas_match(fasta1: str, fasta2: str) -> bool:
    """Check if two fasta files are identical"""
    seqs1 = fasta_to_dict(fasta1)
    seqs2 = fasta_to_dict(fasta2)

    if len(seqs1) != len(seqs2):
        logger.error(f"Number of sequences in {fasta1} and {fasta2} do not match")
        return False

    for name in seqs1:
        if name not in seqs2:
            logger.error(f"Sequence {name} not found in {fasta2}")
            return False
        if seqs1[name] != seqs2[name]:
            logger.error(f"Sequences for {name} do not match")
            return False

    return True


def extract_accession(s: str) -> str:
    m = ACC_REGEX.search(s)
    if m is None:
        logger.error(f"Could not find accession in {s}")
        sys.exit(1)
    try:
        acc = m.group("acc")
    except IndexError:
        logger.error(f"Could not find accession in {s}")
        sys.exit(1)

    return acc


def fetch_checkm_info(
    genomes: List[Path], tmpdir: Path
) -> Dict[str, Tuple[float, float, float]]:
    """Fetch the CheckM info for a list of genomes"""
    accession_list = tmpdir / "accessions.txt"
    accessions = []
    with open(accession_list, "w") as f:
        for genome in genomes:
            p = str(genome)
            acc = extract_accession(p)
            accessions.append(acc)
            print(acc, file=f)

    args = [
        "datasets",
        "summary",
        "genome",
        "accession",
        "--inputfile",
        str(accession_list),
    ]
    logger.debug(f"Running NCBI Datasets with args: {args}")
    proc = subprocess.run(
        args,
        capture_output=True,
        text=True,
    )

    if proc.returncode != 0:
        logger.error("Error running NCBI Datasets")
        logger.error(proc.stderr)
        sys.exit(1)

    data = json.loads(proc.stdout)
    accessions = set(accessions)
    acc_with_report = set()
    reports = data["reports"]
    checkm_info = {}
    for info in reports:
        acc, completeness_percentile, completeness, contamination = (
            extract_checkm_stats(info)
        )
        if acc in accessions:
            acc_with_report.add(acc)
        else:
            logger.warning(f"Accession {acc} not found in the list of genomes")

        checkm_info[acc] = (completeness_percentile, completeness, contamination)

    missing_acc = accessions - acc_with_report
    if missing_acc:
        logger.error(
            f"Could not find reports for the following accessions: {missing_acc}"
        )
        sys.exit(1)

    return checkm_info


def extract_checkm_stats(report: dict) -> Tuple[str, float, float, float]:
    """Extract the CheckM info of a genome from the NCBI assembly summary report"""
    if not report:
        logger.error("Empty report")
        sys.exit(1)

    acc = report.get("accession")
    if acc is None:
        logger.error("Accession not found in report")
        sys.exit(1)

    assembly_status = report.get("assembly_info", {}).get("assembly_status")
    if assembly_status is None:
        logger.error(f"Assembly status not found for {acc}")
        sys.exit(1)
    elif assembly_status.lower() == "suppressed":
        logger.warning(
            f"Assembly status is suppressed for {acc}, assuming 0% completeness and 100% contamination..."
        )
        return acc, 0.0, 0.0, 100.0

    checkm_info = report.get("checkm_info")
    if checkm_info is None:
        logger.warning(
            f"CheckM info not found for {acc}, assuming 100% completeness and 0% contamination..."
        )
        return acc, 100.0, 100.0, 0.0
    try:
        completeness_percentile = float(checkm_info["completeness_percentile"])
    except KeyError:
        logger.error(f"Completeness percentile not found for {acc}")
        sys.exit(1)
    try:
        completeness = float(checkm_info["completeness"])
    except KeyError:
        logger.error(f"Completeness not found for {acc}")
        sys.exit(1)
    try:
        contamination = float(checkm_info["contamination"])
    except KeyError:
        logger.warning(f"Contamination not found for {acc}, using 0.0%")
        contamination = 0.0
    return acc, completeness_percentile, completeness, contamination


def main():
    args = parse_args()
    setup_logging(args.verbose, args.quiet)

    outdir = Path(args.outdir).resolve()
    # Create output directory if it doesn't exist
    if outdir.exists() and args.force:
        logger.info(
            f"Output directory {outdir} already exists, but will be overwritten"
        )
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    dependencies = ["skani", "bcftools", "paftools.js", "datasets"]

    if args.existing is None:
        dependencies.append("genome_updater.sh")

    for dep in dependencies:
        if not is_executable(dep):
            logger.error(f"{dep} is not installed or not executable")
            sys.exit(1)

    check_bcftools_version(MIN_BCFTOOLS_VERSION)

    if args.tmpdir is None:
        tmpdir = Path(tempfile.mkdtemp())
    else:
        tmpdir = Path(args.tmpdir).resolve()
        tmpdir.mkdir(parents=True, exist_ok=True)

    if not args.no_cleanup:
        atexit.register(shutil.rmtree, tmpdir)
    else:
        logger.info(f"Temporary files will be kept in {tmpdir}")

    if args.existing is None:
        logger.info("Downloading genomes...")
        download_genomes(
            str(tmpdir),
            args.database,
            args.taxonomy,
            args.species,
            args.threads,
            args.asmlvl,
            args.max_asm,
        )
        # get all files with extension fna.gz in the tmpdir
        genomes = [p.resolve() for p in tmpdir.rglob("*.fna.gz")]

        if len(genomes) == 0:
            logger.error("No genomes found downloaded")
            sys.exit(1)

        logger.success(f"Downloaded {len(genomes)} genomes")
    else:
        logger.info("Using existing genomes...")
        genomes = []
        # loop over all files in the existing genomes directory
        for p in Path(args.existing).rglob("*"):
            # keep only fasta files
            if FASTA_REGEX.match(p.name):
                genomes.append(p.resolve())

        if len(genomes) == 0:
            logger.error(f"No genomes found in {args.existing}")
            sys.exit(1)

        logger.success(f"Found {len(genomes)} genomes")

    # write the list of genomes to a file
    genomes_files = tmpdir / "genomes.txt"
    genomes_files.write_text("\n".join([str(g) for g in genomes]))

    # fetch CheckM info for each genome
    logger.info("Fetching CheckM info for each genome...")
    checkm_info = fetch_checkm_info(genomes, tmpdir)

    logger.info("Sketching genomes...")
    sketch_dir = tmpdir / "sketches"
    sketch_genomes(str(genomes_files), str(sketch_dir), args.threads)

    logger.info(f"Calculating ANI between the genomes and {args.reference_genome}...")

    dist_mtx = outdir / "ani.tsv"
    calculate_ANI(
        str(sketch_dir), args.reference_genome, str(args.threads), str(dist_mtx)
    )

    # add checkm info to the ANI file
    tmpdist = tmpdir / "ani.tsv.tmp"
    with open(dist_mtx, "r") as f, open(tmpdist, "w") as f_out:
        header = next(f).rstrip().split("\t")
        header.extend(["completeness_percentile", "completeness", "contamination"])
        # remove the Ref_name and Query_name columns
        ref_name_ix = header.index("Ref_name")
        query_name_ix = header.index("Query_name")
        header.pop(ref_name_ix)
        header.pop(query_name_ix - 1)
        print("\t".join(header), file=f_out)
        for line in f:
            fields = line.rstrip().split("\t")
            p = fields[0]
            acc = extract_accession(p)
            fields.pop(ref_name_ix)
            fields.pop(query_name_ix - 1)
            completeness_percentile, completeness, contamination = checkm_info[acc]
            fields.extend(
                [str(completeness_percentile), str(completeness), str(contamination)]
            )
            print("\t".join(fields), file=f_out)

    shutil.move(tmpdist, dist_mtx)

    logger.info(f"You can find the full list of genomes and ANIs in {dist_mtx}")

    with open(dist_mtx, "r") as f:
        distances = []
        header = next(f).rstrip().split("\t")
        ani_ix = header.index("ANI")
        completeness_ix = header.index("completeness")
        contamination_ix = header.index("contamination")
        completeness_percentile_ix = header.index("completeness_percentile")

        for i, line in enumerate(f, start=1):
            fields = line.strip().split("\t")
            ani = float(fields[ani_ix])
            completeness = float(fields[completeness_ix])
            contamination = float(fields[contamination_ix])
            completeness_percentile = float(fields[completeness_percentile_ix])
            acc = extract_accession(fields[0])

            in_ani_range = args.min_ani <= ani <= args.max_ani
            passes_quality_check = (
                completeness >= args.min_completeness
                and contamination <= args.max_contamination
                and completeness_percentile >= args.min_completeness_percentile
            )

            if not in_ani_range:
                logger.trace(
                    f"Genome {acc} has ANI {ani}, which is outside the range {args.min_ani} - {args.max_ani}"
                )
            elif not passes_quality_check:
                logger.trace(
                    f"Genome {acc} has completeness {completeness}%, contamination {contamination}%, and completeness percentile {completeness_percentile}%, which do not pass the quality check"
                )
            else:
                distances.append((ani, i))

    if len(distances) == 0:
        logger.error(
            f"No genomes found that pass ANI and quality criteria: ANI in range {args.min_ani} - {args.max_ani}, completeness >= {args.min_completeness}, contamination <= {args.max_contamination}, and completeness percentile >= {args.min_completeness_percentile}"
        )
        sys.exit(1)

    # get the index of the line closest to ANI
    _, lineno = argclosest(distances, args.ani)

    # get the genome name from the line number
    with open(dist_mtx, "r") as f:
        donor_fields = f.readlines()[lineno].strip().split("\t")
        donor_genome = Path(donor_fields[0]).resolve()
        donor_ani = float(donor_fields[ani_ix])
        donor_completeness = float(donor_fields[completeness_ix])
        donor_contamination = float(donor_fields[contamination_ix])
        donor_completeness_percentile = float(donor_fields[completeness_percentile_ix])

    logger.success(
        f"Using {donor_genome.name}, which has ANI {donor_ani}, {donor_completeness}% completeness, {donor_contamination}% contamination, and {donor_completeness_percentile}% completeness percentile"
    )

    # copy the genome to the output directory and decompress it if necessary
    mutdonor = outdir / "mutdonor.fna"
    with open(mutdonor, "w") as outfile:
        if is_file_gzipped(donor_genome):
            outfile.write(gzip.open(donor_genome, "rb").read().decode())
        else:
            outfile.write(open(donor_genome, "r").read())

    logger.success(f"Mutant donor genome saved to {mutdonor}")

    # do the same for the reference genome
    reference_genome_out = outdir / "reference.fna"
    with open(reference_genome_out, "w") as outfile:
        if is_file_gzipped(args.reference_genome):
            outfile.write(gzip.open(args.reference_genome, "rb").read().decode())
            logger.success(
                f"Reference genome decompressed and saved to {reference_genome_out}"
            )
        else:
            outfile.write(open(args.reference_genome, "r").read())
            logger.success(f"Reference genome saved to {reference_genome_out}")

    reference_genome = reference_genome_out

    logger.info(f"Generating variants between {reference_genome} and {mutdonor}...")

    mm2_vcf = outdir / "minimap2.vcf.gz"
    generate_minimap2_variants(
        str(reference_genome), str(mutdonor), str(mm2_vcf), args.threads
    )

    logger.debug(f"Minimap2 variants saved to {mm2_vcf}")

    dnadiff_vcf = outdir / "dnadiff.vcf.gz"
    generate_dnadiff_variants(
        tmpdir, str(reference_genome), str(mutdonor), str(dnadiff_vcf), args.threads
    )

    logger.debug(f"DNAdiff variants saved to {dnadiff_vcf}")

    logger.info("Merging, normalising, and filtering variants...")

    apply_vcf = outdir / "apply.vcf.gz"
    merge_and_filter_variants(
        tmpdir,
        [str(mm2_vcf), str(dnadiff_vcf)],
        str(reference_genome),
        str(apply_vcf),
        args.max_indel,
        args.remove_overlaps,
    )

    vcfstats_file = outdir / "vcfstats.txt"
    vcfstats(str(apply_vcf), str(vcfstats_file))
    logger.success(f"VCF stats saved to {vcfstats_file}")

    logger.info("Creating mutant reference genome...")

    mutref = outdir / "mutreference.fna"
    create_mutant_reference(
        str(reference_genome),
        str(apply_vcf),
        str(mutref),
    )

    logger.success(f"Mutant reference genome saved to {mutref}")
    logger.success(
        f"VCF used to create the mutant reference genome saved to {apply_vcf}"
    )

    logger.info("Creating truth VCF...")
    inverse_outvcf = tmpdir / "inverse.vcf"
    inverse_vcf(str(apply_vcf), str(inverse_outvcf))

    # compress and index the inverse VCF
    compressed_inverse_vcf = tmpdir / "inverse.vcf.gz"
    cmd = f"bcftools view --write-index -o {compressed_inverse_vcf} {inverse_outvcf}"
    logger.debug(f"Running bcftools view: {cmd}")
    proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    if proc.returncode != 0:
        logger.error(f"bcftools view failed with return code {proc.returncode}")
        logger.error(proc.stderr)
        sys.exit(1)
    else:
        logger.trace(f"bcftools view output:\n{proc.stderr}")

    # check the inverse VCF, applied to the mutant reference, should give the original reference
    inverse_ref = tmpdir / "inverse_reference.fna"
    create_mutant_reference(str(mutref), str(compressed_inverse_vcf), str(inverse_ref))
    seqs_match = fastas_match(str(reference_genome), str(inverse_ref))

    if not seqs_match:
        logger.error("Inverse VCF did not produce the original reference genome")
        sys.exit(1)

    logger.debug("Inverse VCF produced the original reference genome")

    final_vcf = outdir / "truth.vcf.gz"
    shutil.copy(compressed_inverse_vcf, final_vcf)
    # and copy the index
    shutil.copy(str(compressed_inverse_vcf) + ".csi", str(final_vcf) + ".csi")

    logger.success(f"Final truth variants saved to {final_vcf}")
    logger.success("Done!")


if __name__ == "__main__":
    main()
