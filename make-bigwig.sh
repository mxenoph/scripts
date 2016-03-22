#!/usr/bin/env bash

# Adapted from https://github.com/klmr/pol3-seq/blob/master/scripts/bigwig

input="$1"  # BAM file
genome="$2" # Chromosome lenghts for assembly
output="$3" # Bigwig file
outputbg="${output%%.bw}.bedgraph"  # Same as "${output/.bw/.bedgraph}"

bedtools genomecov -bg -ibam "$input" -g "$genome" > "$outputbg"
bedGraphToBigWig "$outputbg" "$genome" "$output"
