#!/bin/env bash

while read line
do
    echo $line >> $2
    bigWigCorrelate $line >> $2
done < $1
