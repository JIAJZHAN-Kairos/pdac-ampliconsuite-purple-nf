#!/usr/bin/env bash
# Check that a BAM is readable, indexed, and uses the GRCh38 primary contig lengths.
set -euo pipefail
bam=$1
fai=$2
samtools quickcheck -v "$bam"
test -s "$bam.bai"
samtools idxstats "$bam" > "$bam.index.tsv"
awk 'NR==FNR { if ($1 ~ /^chr([0-9]+|X|Y)$/) ref[$1]=$2; next }
     $1 in ref { if ($2 != ref[$1]) exit 1; delete ref[$1] }
     END { if (length(ref)) exit 1 }' "$fai" "$bam.index.tsv"
