#!/usr/bin/env bash

input="$1"
output="$2"

samtools view -h -o $output $input
