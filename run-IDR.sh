#! /usr/bin/env bash

# For combining the files for IDR comparison
# e.g. for rep1 rep2 rep3 it will print "rep1 rep2"
# "rep1 rep3" and "rep2 rep3"

#run-IDR.sh #{{{
usage() { echo "Usage: $0
                [-p <relaxed macs2 peaks (non truncated list)>] 
                [-d <relaxed peaks from pooled>]
                [-s <peaks on pseudoreplicates from pooled experiment> ]
                [-e <self-pseudoreplicates peaks> ]
                [-o <output directory> ]
                [-n <options> ]
                [-t <top number of overlapping peaks to select> ]
                [-a <top number of peaks to keep for pseudoreplicates> ]
                [-b <top number of peaks to keep for pooled pseudoreplicates> ]
                [-q <IDR threshold>]
                [-f <IDR threshold for pseudoreplicates> ]
                [-g <IDR threshold for pooled pseudoreplicates> ]" 1>&2; exit 1; }

# : means takes an argument but not mandatory (if mandatory will have to check after)
options=':p:d:s:e:o:n:t:a:b:q:f:g:h'
while getopts $options option
do
    case $option in
        p ) peaks=(${OPTARG}) ;;
        d ) pooled=(${OPTARG}) ;;
        s ) pooled_ps=(${OPTARG}) ;;
        e ) pseudoreps=(${OPTARG}) ;;
        o ) output_path=${OPTARG} ;;
        n ) OPTIONS=${OPTARG} ;;
        t ) TOP=${OPTARG} ;;
        a ) ps_TOP=${OPTARG} ;;
        b ) pps_TOP=${OPTARG} ;;
        q ) IDR=${OPTARG} ;;
        f ) ps_IDR=${OPTARG} ;;
        g ) pps_IDR=${OPTARG} ;;
        h ) usage ;;
        : ) echo "Missing option argument for -$OPTARG" >&2; usage ;;
        \?) echo "Unknown option: -$OPTARG" >&2; usage ;;
        * ) usage ;;
    esac
done
shift $((OPTIND-1))

if [ -z "$OPTIONS" ]
then
    OPTIONS='0 F p.value'
fi

# if the number of for truncating the macs2 relaxed peaks is not provided
# use top 100K by default (underlying assumption: it's a big genome)
if [ -z "$TOP" ]
then
    TOP=150000
fi

# if IDR threshold not provided use 0.02 by default
if [ -z "$IDR" ]
then
    IDR=0.02
fi

#}}}

# Path should be the same for all files used in comparisons#{{{
# It definetely is if data from makefile
mkdir -p ${output_path}/batch-consistency
mkdir -p ${output_path}/thresholds
mkdir -p ${output_path}/final-set

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
        # the wrong order.
        sort -k 8nr,8nr ${peaks[$((i-1))]} | head -n ${TOP} >  ${peaks[$((i-1))]/_peaks.narrowPeak/_top$(($TOP/1000))K}
        sort -k 8nr,8nr ${peaks[$((j-1))]} | head -n ${TOP} >  ${peaks[$((j-1))]/_peaks.narrowPeak/_top$(($TOP/1000))K}
        echo  ${peaks[$((i-1))]/_peaks.narrowPeak/_top$(($TOP/1000))K}
        
        base_one=$(basename ${peaks[$((i-1))]} | sed 's/.*-//g' | sed 's/_peaks.*//g')
        base_two=$(basename ${peaks[$((j-1))]} | sed 's/.*-//g' | sed 's/_peaks.*//g')
        current_comparison=${output_path}/batch-consistency/${target}-${base_one}vs${base_two}-${IDR}
        echo "ORIGINAL REPS: ${current_comparison}"

#        Rscript ~/local/idrCode/batch-consistency-analysis.r ${peaks[$((i-1))]/_peaks.narrowPeak/_top$(($TOP/1000))K} ${peaks[$((j-1))]/_peaks.narrowPeak/_top$(($TOP/1000))K} -1 "${current_comparison}" ${OPTIONS}
        comparisons+=("${current_comparison}")
        tmp=$(awk -v idr=${IDR} '$11 <= idr {print $0}' "${current_comparison}-overlapped-peaks.txt" | wc -l)
        num_peaks+=(${tmp})

       if [ -e ${output_path}/thresholds/${target}-${base_one/_*}_original-replicate-thresholds.tsv ] && [ "$i" -eq 1 ]
       then
           # Remove file because thresholds are appended to it
           echo "Removing ${output_path}/thresholds/${target}-${base_one/_*}_original-replicate-thresholds.tsv"
           rm ${output_path}/thresholds/${target}-${base_one/_*}_original-replicate-thresholds.tsv
       fi
#        echo -e "${current_comparison}\t${tmp}" >> ${output_path}/thresholds/${target}-${base_one/_*}_original-replicate-thresholds.tsv
        num_peaks+=($(awk -v idr=${IDR} '$11 <= idr {print $0}' "${output_path}/batch-consistency/${target}-${base_one}vs${base_two}-${IDR}-overlapped-peaks.txt" | wc -l))
    done
done

#Rscript ~/local/idrCode/batch-consistency-plot.r ${#comparisons[@]} "${output_path}/batch-consistency/${target}-${base_one/_*}_reps" ${comparisons[*]}

IFS=$'\n'
# Original replicate threshold
original_replicate_thr=$(echo "${num_peaks[*]}" | sort -nr | head -n1)
echo "ORIGINAL REPLICATE THRESHOLD: ${original_replicate_thr}"

#sort -k 8nr,8nr ${pooled} | head -n ${original_replicate_thr} > ${output_path}/final-set/$(sed "s/_[^-]*$/_conservative_${IDR}.narrowPeak/g" <<< $(basename ${pooled}))
#}}}

# On pseudoreplicates #{{{
if [ ! -z "$pseudoreps" ] || [ ! -z "$pooled_ps" ]
then
    echo 'IDR on pseudoreplicates'
    if [ -z "$ps_TOP" ]
    then
        ps_TOP=150000
    fi

    if [ -z "$ps_IDR" ]
    then
        ps_IDR=0.02
    fi

    if [ -z "$pps_TOP" ]
    then
        ps_TOP=150000
    fi

    if [ -z "$pps_IDR" ]
    then
        ps_IDR=0.02
    fi

    IFS=' ' read -r -a ps00 <<< ${pseudoreps[@]//*pseudoreps01*/}
    IFS=' ' read -r -a ps01 <<< ${pseudoreps[@]//*pseudoreps00*/}

    target=$(basename ${ps00[0]} | sed 's/-[^-]*$//g')
    comparisons=()
    num_pseudoreps=()

    # just constructing the output name for the self-consistency thresholds file
    output_thr=$(basename ${ps00[0]} | sed 's/.*-//g' | sed 's/_.*//g')
    output_thr="${output_path}/thresholds/${target}-${output_thr}_self-consistency-thresholds.tsv"

    if [ -e ${output_thr} ]
    then
        # Remove file because thresholds are appended to it
        echo "Removing ${output_thr}"
        rm ${output_thr}
    fi

    for i in `seq 1 ${#ps00[@]}`
    do
        # truncate list of relaxed pseudoreps
        # ToDO: this is hard-coded so if p-value not used as threshold for the peak calling it will be 
        # the wrong order. Either 
#        sort -k 8nr,8nr ${ps00[$((i-1))]} | head -n ${ps_TOP} >  ${ps00[$((i-1))]/_peaks.narrowPeak/_top$(($ps_TOP/1000))K}
#        sort -k 8nr,8nr ${ps01[$((i-1))]} | head -n ${ps_TOP} >  ${ps01[$((i-1))]/_peaks.narrowPeak/_top$(($ps_TOP/1000))K}
            
        base_one=$(basename ${ps00[$((i-1))]} | sed 's/.*-//g' | sed 's/_peaks.*//g')
        base_two=$(basename ${ps01[$((i-1))]} | sed 's/.*-//g' | sed 's/_peaks.*//g')
        current_comparison=${output_path}/batch-consistency/${target}-${base_one}vs${base_two}-${ps_IDR}
        echo ${current_comparison}

        Rscript ~/local/idrCode/batch-consistency-analysis.r ${ps00[$((i-1))]/_peaks.narrowPeak/_top$(($ps_TOP/1000))K} ${ps01[$((i-1))]/_peaks.narrowPeak/_top$(($ps_TOP/1000))K} -1 "${current_comparison}" ${OPTIONS}
        ls ${current_comparison}*
        comparisons+=("${current_comparison}")
        tmp=$(awk -v idr=${ps_IDR} '$11 <= idr {print $0}' "${current_comparison}-overlapped-peaks.txt" | wc -l)
        num_pseudoreps+=(${tmp})
        echo -e "${current_comparison}\t${tmp}" >> ${output_path}/thresholds/${target}-${base_one/_*}_self-consistency-thresholds.tsv

    done
    exit
    Rscript ~/local/idrCode/batch-consistency-plot.r ${#comparisons[@]} "${output_path}/batch-consistency/${target}-${base_one/_*}_pseudoreps" ${comparisons[*]}

    echo "IDR on pooled pseudoreplicates"
    sort -k 8nr,8nr ${pooled_ps[0]} | head -n ${pps_TOP} >  ${pooled_ps[0]/_peaks.narrowPeak/_top$(($pps_TOP/1000))K}
    sort -k 8nr,8nr ${pooled_ps[1]} | head -n ${pps_TOP} >  ${pooled_ps[1]/_peaks.narrowPeak/_top$(($pps_TOP/1000))K}
    
    base_one=$(basename ${pooled_ps[0]} | sed 's/.*-//g' | sed 's/_peaks.*//g')
    base_two=$(basename ${pooled_ps[1]} | sed 's/.*-//g' | sed 's/_peaks.*//g')
    current_comparison=${output_path}/batch-consistency/${target}-${base_one}vs${base_two}-${pps_IDR}
    Rscript ~/local/idrCode/batch-consistency-analysis.r ${pooled_ps[0]/_peaks.narrowPeak/_top$(($pps_TOP/1000))K} ${pooled_ps[1]/_peaks.narrowPeak/_top$(($pps_TOP/1000))K} -1 "${current_comparison}" ${OPTIONS}
    Rscript ~/local/idrCode/batch-consistency-plot.r 1 "${output_path}/batch-consistency/${target}-${base_one/_*}_pooled_pseudoreps" ${current_comparison}

    pooled_replicate_thr=$(awk -v idr=${pps_IDR} '$11 <= idr {print $0}' "${current_comparison}-overlapped-peaks.txt" | wc -l)
    
    # just constructing the output name for the pooled-consistency thresholds file
    if [ -e "${output_path}/thresholds/${target}-${base_one/_*}_pooled-consistency-thresholds.tsv" ]
    then
        # Remove file because thresholds are appended to it
        echo "Removing ${output_path}/thresholds/${target}-${base_one/_*}_pooled-consistency-thresholds.tsv"
        rm ${output_path}/thresholds/${target}-${base_one/_*}_pooled-consistency-thresholds.tsv
    fi

    echo -e "${current_comparison}\t${pooled_replicate_thr}" >> ${output_path}/thresholds/${target}-${base_one/_*}_pooled-consistency-thresholds.tsv
    echo "POOLED REPLICATE THRESHOLD: ${pooled_replicate_thr}"

    IFS=$'\n'
    # Pooled pseudoreplicate threshold
    optimal_thr=$(cut -f 2 ${output_path}/thresholds/${target}-${base_one/_*}_pooled-consistency-thresholds.tsv ${output_path}/thresholds/${target}-${base_one/_*}_original-replicate-thresholds.tsv | sort -nr | head -n 1)
    echo "OPTIMAL THRESHOLD: ${optimal_thr}"

    sort -k 8nr,8nr ${pooled} | head -n ${optimal_thr} > ${output_path}/final-set/$(sed "s/_[^-]*$/_optimal_${IDR}.narrowPeak/g" <<< $(basename ${pooled}))

fi
#}}}
