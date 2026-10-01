"""Clip reads to a reference window, from `samtools view` text on stdin.

    samtools view -F 0x904 aln.bam contig:START-END | clip_reads.py START END > clipped.fastq

Each primary alignment overlapping the 1-based inclusive window [START, END] is cut down to
the bases that align inside it, then put back in the read's sequenced orientation. Every
position in the window is then covered by the same reads it would have in the whole genome
(no ramp at the window's ends, as there is when only reads lying wholly inside are kept),
so depth-capped subsampling of the fixture behaves as it does on a full Read set.
"""

import re
import sys

start, end = int(sys.argv[1]) - 1, int(sys.argv[2])  # 0-based, half-open
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
CIGAR_OP = re.compile(r"(\d+)([MIDNSHP=X])")

for line in sys.stdin:
    f = line.rstrip("\n").split("\t")
    name, flag, pos, cigar, seq, qual = f[0], int(f[1]), int(f[3]) - 1, f[5], f[9], f[10]
    r, q, lo, hi = pos, 0, None, None
    for n, op in CIGAR_OP.findall(cigar):
        n = int(n)
        if op in "M=X":
            a, b = max(r, start), min(r + n, end)
            if a < b:
                lo = q + (a - r) if lo is None else lo
                hi = q + (b - r)
            r, q = r + n, q + n
        elif op in "DN":
            r += n
        elif op in "IS":
            q += n
    if lo is None:
        continue
    s, ql = seq[lo:hi], qual[lo:hi]
    if flag & 16:
        s, ql = s.translate(COMP)[::-1], ql[::-1]
    sys.stdout.write(f"@{name}\n{s}\n+\n{ql}\n")
