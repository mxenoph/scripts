#!/usr/bin/env bash

target_dir="$(dirname "$1")"/fastqc
target_base="$(basename "$1")"
mkdir -p "$target_dir"

/homes/mxenoph/local/fastqc -o "$target_dir" $1
