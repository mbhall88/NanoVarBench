#!/usr/bin/env bash
set -euo pipefail

# This script will produce a list of genomes within a given mash distance of the provided reference genome for the given species
# Function to test if a dependency is available
check_dependency() {
    command -v "$1" >/dev/null 2>&1 || {
        echo >&2 "Dependency '$1' is required but not installed. Aborting."
        exit 1
    }
}

# function to make sure bcftools version is >=1.18
check_bcftools_version() {
    bcftools_version=$(bcftools --version | grep -oP "(?<=bcftools )\d+\.\d+")
    if [ "$(echo "$bcftools_version >= 1.18" | bc)" -eq 0 ]; then
        echo "bcftools version must be >=1.18" >&2
        exit 1
    fi
}

db="refseq"
organism_grp="bacteria"
file_type="genomic.fna.gz"
taxonomy="gtdb"
threads=1
sketch_size=10000
asmlvl="complete genome"
outdir=$(realpath "mutref")
tmpdir=""
species=""
reference_genome=""
existing_genomes=""
# Set the default mash distance to 0.005
mash_distance=0.005
max_indel=50
cleanup_tmp_files=true

# Usage statement
usage() {
    cat <<EOF
Usage: $0 -r <reference_genome> -s <species> [-d <mash_distance>] [-c] [-t <threads>] [-l <assembly_level>] [-T <taxonomy>] [-D <database>] [-S <sketch_size>] [-o <outdir>] [-e <tmpdir>] [-f <existing_genomes>] [-I <max_indel] [-h]
    -r: reference genome
    -s: species - e.g., "s__Escherichia coli" (GTDB) or "Escherichia coli" (NCBI)
    -d: pick genome with mash distance closest to this value (default: $mash_distance)
    -c: do not cleanup temporary files
    -t: number of threads (default: $threads)
    -l: assembly level. (comma-separated entries, empty for all) [complete genome, chromosome, scaffold, contig] (default: $asmlvl)
    -T: taxonomy database. [gtdb, ncbi] (default: $taxonomy)
    -D: database. (comma-separated entries) [refseq, genbank] (default: $db)
    -S: mash sketch size (default: $sketch_size)
    -o: output directory (default: mutref)
    -e: temporary directory
    -f: directory containing existing genomes (skips downloading genomes)
    -I: maximum indel size (default: $max_indel)
    -h: show this message
EOF
}

# Define the command line options
while getopts "r:s:d:ct:l:T:D:S:o:e:h" opt; do
    case $opt in
    r)
        reference_genome="$OPTARG"
        ;;
    s)
        species="$OPTARG"
        ;;
    d)
        mash_distance="$OPTARG"
        ;;
    c)
        cleanup_tmp_files=false
        ;;
    t)
        threads="$OPTARG"
        ;;
    l)
        asmlvl="$OPTARG"
        ;;
    T)
        taxonomy="$OPTARG"
        ;;
    D)
        db="$OPTARG"
        ;;
    S)
        sketch_size="$OPTARG"
        ;;
    o)
        outdir=$(realpath "$OPTARG")
        ;;
    e)
        tmpdir="$OPTARG"
        ;;
    f)
        existing_genomes="$OPTARG"
        ;;
    h)
        usage
        exit 0
        ;;
    \?)
        echo "Invalid option -$OPTARG" >&2
        usage
        exit 1
        ;;
    esac
done

# Check if required options are set
if [ -z "$reference_genome" ]; then
    echo "Error: reference genome (-r) not provided" >&2
    usage
    exit 1
fi

if [ -z "$species" ] && [ -z "$existing_genomes" ]; then
    echo "Error: species (-s) not provided" >&2
    usage
    exit 1
fi

# List of dependencies
dependencies=("mash" "fd" "python" "varifier" "bcftools" "paftools.js" "gzip")

# add genome_updater.sh to the list of dependencies if existing_genomes is not set
if [ -z "$existing_genomes" ]; then
    dependencies+=("genome_updater.sh")
fi

# Check if each dependency is available
for dependency in "${dependencies[@]}"; do
    check_dependency "$dependency"
done

check_bcftools_version

# if tmpdir is not set, create a temporary directory
if [ -z "$tmpdir" ]; then
    tmpdir=$(mktemp -d)
fi

if $cleanup_tmp_files; then
    trap 'rm -rf "$tmpdir"' EXIT
else
    echo "Temporary files will be kept in $tmpdir" >&2
fi

if [ ! -d "$outdir" ]; then
    mkdir -p "$outdir"
fi

# if existing_genomes is not set, download the genomes
if [ -z "$existing_genomes" ]; then
    echo "Downloading the species genomes..." >&2
    genome_updater.sh -d "$db" -g $organism_grp -f $file_type -o "$tmpdir" -M "$taxonomy" -T "$species" -m -a -t "$threads" -l "$asmlvl"
    fd -a -t f -e fna.gz . "$tmpdir"
else
    echo "Using existing genomes in $existing_genomes" >&2
    fd -a -t f -e fna.gz . "$existing_genomes"
fi >"$tmpdir"/genomes.txt

if [ -z "$existing_genomes" ]; then
    echo "Downloaded $(wc -l <"$tmpdir"/genomes.txt) genomes" >&2
else
    echo "Found $(wc -l <"$tmpdir"/genomes.txt) existing genomes" >&2
fi

echo "Finding genomes around a mash distance of $mash_distance..." >&2

mash sketch -p "$threads" -o "$tmpdir"/sketch -l "$tmpdir"/genomes.txt -s "$sketch_size"

dist_mtx="$outdir"/distances.tsv
mash dist -p "$threads" -s "$sketch_size" "$tmpdir"/sketch.msh "$reference_genome" | sort -g -k3 >"$dist_mtx"

# get the index of the line closest to mash_distance
lineno=$(cut -f 3 "$dist_mtx" | python -c "import sys;dists=[float(i) for i in sys.stdin.readlines()];print(min(enumerate(dists), key=lambda t:abs(t[1]-0.005))[0]+1)")

# get the genome name from the line number
genome=$(sed -n "${lineno}p" "$dist_mtx" | cut -f 1)
dist=$(sed -n "${lineno}p" "$dist_mtx" | cut -f 3)
shared_hashes=$(sed -n "${lineno}p" "$dist_mtx" | cut -f 5)

echo "Using $(basename "$genome"), which has a mash distance of $dist and $shared_hashes shared hashes" >&2

# copy the genome to the output directory and decompress it if necessary
mutdonor="$outdir"/mutdonor.fna
if file -b "$genome" | grep -q "gzip compressed data"; then
    gzip -dc "$genome"
else
    cat "$genome"
fi >"$mutdonor"

echo "Mutant donor genome saved to $mutdonor" >&2

# do the same for the reference genome
if file -b "$reference_genome" | grep -q "gzip compressed data"; then
    gzip -dc "$reference_genome"
else
    cat "$reference_genome"
fi >"$outdir"/reference.fna
reference_genome="$outdir"/reference.fna

echo "Generating variants between $reference_genome and $mutdonor..." >&2
varifier_out="$tmpdir"/varifier
varifier make_truth_vcf "$mutdonor" "$reference_genome" "$varifier_out" --cpus "$threads"

vcf="$outdir"/ref_donor.vcf.gz
bcftools view -e "abs(ILEN)>$max_indel" -o "$vcf" --write-index "$varifier_out"/04.truth.vcf
paftools.js vcfstat "$vcf" | tee "$outdir"/ref_donor.vcfstat
echo "Truth VCF saved to $vcf" >&2

echo "Applying variants to $reference_genome..." >&2
bcftools consensus -f "$reference_genome" "$vcf" >"$outdir"/mutreference.fna
echo "Mutant reference genome saved to $outdir/mutreference.fna" >&2

echo "Done!" >&2