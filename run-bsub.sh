#!/usr/bin/env bash

# Adapted from https://github.com/klmr/pol3-seq/blob/master/scripts/bsub

args=("$@")
# Last argument provided is always the job I want to run
job=("${args[${#args[@]}-1]}")
# Delete the job from the array of bsub arguments
unset args[${#args[@]}-1]
# Job name is anything before the first space thus the application running
jobname="$(basename "$( sed 's/ .*$//' <<< "$job" )")"
# what if logs dir not present? it's not created in the makefile
# How do you get the job id?
logfile="logs/$jobname-%J.log"
errfile="logs/$jobname-%J.err"

bsub -J "$jobname" -o "$logfile" -e "$errfile" ${args[@]} "$job"

