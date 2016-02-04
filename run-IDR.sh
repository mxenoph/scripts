#! /usr/bin/env bash

# For combining the files for IDR comparison
# e.g. for rep1 rep2 rep3 it will print "rep1 rep2"
# "rep1 rep3" and "rep2 rep3"

#run-IDR.sh #{{{
usage() { echo "Usage: $0 [-p <relaxed macs2 peaks (non truncated list)>] [-o <output directory>] [-s <options>] [-t <top number of overlapping peaks to select>] [-q <IDR threshold>]" 1>&2; exit 1; }

# : means takes an argument but not mandatory (if mandatory will have to check after)
options=':p:d:o:s:t:q:h'
while getopts $options option
do
    case $option in
        p ) peaks=(${OPTARG}) ;;
        d ) pooled=(${OPTARG}) ;;
        o ) output_path=${OPTARG} ;;
        s ) OPTIONS=${OPTARG} ;;
        t ) TOP=${OPTARG} ;;
        q ) IDR=${OPTARG} ;;
        h ) usage ;;
        : ) echo "Missing option argument for -$OPTARG" >&2; usage ;;
        \?) echo "Unknown option: -$OPTARG" >&2; usage ;;
        * ) usage ;;
    esac
done
shift $((OPTIND-1))

if [ -z "$peaks" ] || [ -z "$output_path" ] || [ -z "$pooled" ]
then
    usage
fi

if [ -z "$OPTIONS" ]
then
    OPTIONS='0 F p.value'
fi

# if the number of for truncating the macs2 relaxed peaks is not provided
# use top 100K by default (underlying assumption: it's a big genome)
if [ -z "$TOP" ]
then
    TOP=100000
fi

# if IDR threshold not provided use 0.02 by default
if [ -z "$IDR" ]
then
    IDR=0.02
fi
#}}}

# Path should be the same for all files used in comparisons
# It definetely is if data from makefile
target=$(basename ${peaks[0]} | sed 's/-[^-]*$//g')
comparisons=()
num_peaks=()
for i in `seq 1 $((${#peaks[@]}-1))`
do
    x=$((i+1));
    for j in `seq $x ${#peaks[@]}`
    do
        # truncate list of relaxed peaks
        # ToDO: this is hard-coded so if p-value not used as threshold for the peak calling it will be 
        # the wrong order. Either 
        sort -k 8nr,8nr ${peaks[$((i-1))]} | head -n ${TOP} >  ${peaks[$((i-1))]/_peaks.narrowPeak/_top$(($TOP/1000))K}
        sort -k 8nr,8nr ${peaks[$((j-1))]} | head -n ${TOP} >  ${peaks[$((j-1))]/_peaks.narrowPeak/_top$(($TOP/1000))K}
        
        base_one=$(basename ${peaks[$((i-1))]} | sed 's/.*-//g' | sed 's/_peaks.*//g')
        base_two=$(basename ${peaks[$((j-1))]} | sed 's/.*-//g' | sed 's/_peaks.*//g')

        Rscript ~/local/idrCode/batch-consistency-analysis.r ${peaks[$((i-1))]/_peaks.narrowPeak/_top$(($TOP/1000))K} ${peaks[$((j-1))]/_peaks.narrowPeak/_top$(($TOP/1000))K} -1 "${output_path}/${target}-${base_one}vs${base_two}-${IDR}" ${OPTIONS}
        comparisons+=("${output_path}/${target}-${base_one}vs${base_two}-${IDR}")
        num_peaks+=($(awk -v idr=${IDR} '$11 <= idr {print $0}' "${output_path}/${target}-${base_one}vs${base_two}-${IDR}-overlapped-peaks.txt" | wc -l))
    done
done

#echo "Batch consistency analysis on ${path}" > "${output_path}/${target}-${base_one/_*}.idr"

Rscript ~/local/idrCode/batch-consistency-plot.r ${#comparisons[@]} "${output_path}/${target}-${base_one/_*}_reps" ${comparisons[*]}

IFS=$'\n'
max_num_peaks=$(echo "${num_peaks[*]}" | sort -nr | head -n1)

sort -k 8nr,8nr ${pooled} | head -n ${max_num_peaks} > ${output_path}/$(sed "s/_[^-]*$/_conservative_${IDR}.narrowPeak/g" <<< $(basename ${pooled}))
