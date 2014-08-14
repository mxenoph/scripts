#!/usr/bin/env bash

input=("$@")

target_dir="$(dirname "${input[1]}")"
target_base="$(basename "${input[1]}")"
output="$target_dir/$target_base.bam"

samtools merge -f "${output}" "${input}"



