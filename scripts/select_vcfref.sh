#!/usr/bin/env bash
set -euo pipefail

# This script will produce a list of genomes within a given mash distance of the provided reference genome for the given species

db="refseq"
organism_grp="bacteria"
file_type="genomic.fna.gz"
taxonomy="gtdb"
threads=1
sketch_size=10000
asmlvl="complete genome"
output="-"
tmpdir=""
species=""
reference_genome=""
# Set the default mash distance to 0.005
mash_distance=0.005
cleanup_tmp_files=true

# Usage statement
usage() {
    cat <<EOF
Usage: $0 -r <reference_genome> -s <species> [-d <mash_distance>] [-c] [-t <threads>] [-l <assembly_level>] [-T <taxonomy>] [-D <database>] [-S <sketch_size>] [-o <output>]
    -r: reference genome
    -s: species - e.g., "Escherichia coli" or "s__Escherichia coli"
    -d: pick genome with mash distance closest to this value (default: $mash_distance)
    -c: do not cleanup temporary files
    -t: number of threads (default: $threads)
    -l: assembly level. (comma-separated entries, empty for all) [complete genome, chromosome, scaffold, contig] (default: $asmlvl)
    -T: taxonomy database. [gtdb, ncbi] (default: $taxonomy)
    -D: database. (comma-separated entries) [refseq, genbank] (default: $db)
    -S: mash sketch size (default: $sketch_size)
    -o: output file (default: stdout)
    -e: temporary directory
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
        output="$OPTARG"
        ;;
    e)
        tmpdir="$OPTARG"
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
if [ -z "$reference_genome" ] || [ -z "$species" ]; then
    echo "Error: required options -r and -s are not set" >&2
    usage
    exit 1
fi

# Function to test if a dependency is available
check_dependency() {
    command -v "$1" >/dev/null 2>&1 || {
        echo >&2 "Dependency '$1' is required but not installed. Aborting."
        exit 1
    }
}

# List of dependencies
dependencies=("genome_updater.sh" "mash" "fd" "python")

# Check if each dependency is available
for dependency in "${dependencies[@]}"; do
    check_dependency "$dependency"
done

# if tmpdir is not set, create a temporary directory
if [ -z "$tmpdir" ]; then
    tmpdir=$(mktemp -d)
fi

if $cleanup_tmp_files; then
    trap 'rm -rf "$tmpdir"' EXIT
else
    echo "Temporary files will be kept in $tmpdir" >&2
fi

echo "Downloading the species genomes..." >&2

genome_updater.sh -d "$db" -g $organism_grp -f $file_type -o "$tmpdir" -M "$taxonomy" -T "$species" -m -a -t "$threads" -l "$asmlvl"

fd -a -t f -e fna.gz . "$tmpdir" >"$tmpdir"/genomes.txt

echo "Downloaded $(wc -l <"$tmpdir"/genomes.txt) genomes" >&2

echo "Finding genomes around a mash distance of $mash_distance..." >&2

mash sketch -p "$threads" -o "$tmpdir"/sketch -l "$tmpdir"/genomes.txt -s "$sketch_size"

mash dist -p "$threads" -s "$sketch_size" "$tmpdir"/sketch.msh "$reference_genome" | sort -g -k3 >"$tmpdir"/distances.txt

# get the index of the line closest to mash_distance
lineno=$(cut -f 3 "$tmpdir"/distances.txt | python -c "import sys;dists=[float(i) for i in sys.stdin.readlines()];print(min(enumerate(dists), key=lambda t:abs(t[1]-0.005))[0]+1)")

# get the genome name from the line number
genome=$(sed -n "${lineno}p" "$tmpdir"/distances.txt | cut -f 1)
dist=$(sed -n "${lineno}p" "$tmpdir"/distances.txt | cut -f 3)
shared_hashes=$(sed -n "${lineno}p" "$tmpdir"/distances.txt | cut -f 5)

echo "Using $(basename "$genome"), which has a mash distance of $dist and $shared_hashes shared hashes" >&2

if [ "$output" == "-" ]; then
    cat "$genome"
else
    cp "$genome" "$output"
    echo "Copied $genome to $output" >&2
fi
