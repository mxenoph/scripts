#!/usr/bin/env bash

input=$1
output=${input%.*}.insert_sizes

samtools view -f66 $1 | cut -f 9 > ${output}
