#!/usr/bin/env bash

target_dir="$(dirname "$1")"/fastqc
target_base="$(basename "$1")"
mkdir -p "$target_dir"

fastqc -o "$target_base" $1
