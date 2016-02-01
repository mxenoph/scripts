#!/usr/bin/env bash

BROAD=false;
RELAXED=false;
#run-macs2.sh #{{{
usage() { echo "Usage: $0 [-i <IP bam>] [-c <Input bam>] [-o <output_path>] [-g <genome e.g.mm>] [-b <broad-boolean>] [-r <relaxed thresholds for IDR method>]" 1>&2; exit 1; }

# : means takes an argument but not mandatory (if mandatory will have to check after)
options=':i:c:g:o:a:rbh'
while getopts $options option
do
    case $option in
        i ) ip=(${OPTARG}) ;;
        c ) CTRL=${OPTARG} ;;
        g ) genome=${OPTARG} ;;
        o ) OUT=${OPTARG} ;;
        b ) BROAD=true ;;
        a ) ASSEMBLY=${OPTARG} ;;
        r ) RELAXED=true ;;
        h ) usage ;;
        : ) echo "Missing option argument for -$OPTARG" >&2; usage ;;
        \?) echo "Unknown option: -$OPTARG" >&2; usage ;;
        * ) usage ;;
    esac
done
shift $((OPTIND-1))

# Check that mandatory arguments were provided
if [ -z "$ip" ] || [ -z "$genome" ]
then
    usage
fi
#}}}

base="$(basename "$ip")"
# Removing file extension and add macs2
target="${base%.*}"
#if [ ${RELAXED} == true ]
#then
#    target="${base%.*}_IDR"
#fi

options=''

# make assembly mandatory and infer genome size/organism initials from that
#genome=$(sed 's/[0-9]//g' <<< ${assembly} | tr '[:upper:]' '[:lower:]')


# If control provided then will pass it to macs#{{{
if [ ! -z "$CTRL" ]
then
    # --bdg will output also bedgraph files for ip and control and --SPMR specifies that the
    # fragment profiles for the pileup are per million reads
    options="$options -c $CTRL --bdg --SPMR"
    if [ ${RELAXED} == true ]
    then
        options="$options --to-large"
    fi

fi
#}}}

# Controlling output directory#{{{
if [ ! -z "$OUT" ]
then
    if [[ ! $OUT =~ \/macs2$ ]]
    then 
        OUT="$(readlink -m ${OUT})/macs2"
        #mkdir -p $OUT
    fi
    #options="$options --outdir $OUT"
else
    OUT=$PWD/macs2
    #options="$options --outdir $PWD"
fi
#}}}

# Controlling peak shape#{{{
if [ ${BROAD} == true ]
then
    OUT="$OUT/broad"
    if [ ${RELAXED} == false ]
    then
        options="$options --broad -q 0.05 --broad-cutoff 0.05"
    else
        # This uses p-values as the cutoff fro calling peaks. Should only use this option in
        # an IDR analysis.
        options="$options --broad -p 0.001 --broad-cutoff 0.005"
    fi
else
    OUT="$OUT/sharp"
    if [ ${RELAXED} == false ]
    then
        options="$options -q 0.01 --call-summits"
    else
        # This uses p-values as the cutoff fro calling peaks. Should only use this option in
        # an IDR analysis.
        options="$options -p 0.001"
    fi

fi
#}}}

#mkdir -p $OUT
options="$options --outdir $OUT"
# right version of macs is in the virtualenv
#Run macs
macs2 callpeak \
    -t ${ip[*]} \
    --g=$genome -n=$target \
    --verbose 0 $options

if [ ${RELAXED} == true ]
then
    echo 'Skipping pileup'
elif [[ !( -z "$CTRL" ) &&  !( -z "$ASSEMBLY" ) ]]
then
    ASSEMBLY=$(echo $ASSEMBLY | tr '[:upper:]' '[:lower:]')
    macs2 bdgcmp -t "${OUT}/${target}_treat_pileup.bdg" -c "${OUT}/${target}_control_lambda.bdg" \
        -m subtract --outdir $OUT --o-prefix ${target}
    macs2 bdgcmp -t "${OUT}/${target}_treat_pileup.bdg" -c "${OUT}/${target}_control_lambda.bdg" \
        -m logFE --outdir $OUT --o-prefix ${target} -p 0.00001

    chrom_sizes=$(grep ${ASSEMBLY} /nfs/research2/bertone/user/mxenoph/genome_dir/assemblies-annotations.config | cut -f 3)
    bedGraphToBigWig ${OUT}/${target}_subtract.bdg ${chrom_sizes} ${OUT}/${target}_subtract.bw
    bedGraphToBigWig ${OUT}/${target}_logFE.bdg ${chrom_sizes} ${OUT}/${target}_logFE.bw
    
    # Adding track line to broadPeak
    if [ ${BROAD} == true ]
    then
        sed -i "1s/^/track type=broadPeak visibility=3 db=${ASSEMBLY} name=\"${target}_peaks.broadPeak\" description=\"${target} ${ASSEMBLY}\"\n/" ${OUT}/${target}_peaks.broadPeak
        #Rscript ~/source/convert-bed-plus-to-gtf.R -b ${OUT}/${target}_peaks.broadPeak -a ${ASSEMBLY} &
    fi

    #Rscript ~/source/convert-bed-plus-to-gtf.R -b ${OUT}/${target}_peaks.narrowPeak -a ${ASSEMBLY} &

    # Do not remove these files as 
    bedGraphToBigWig ${OUT}/${target}_treat_pileup.bdg ${chrom_sizes} ${OUT}/${target}_treat_pileup.bw
    bedGraphToBigWig ${OUT}/${target}_control_lambda.bdg ${chrom_sizes} ${OUT}/${target}_control_lambda.bw
    echo rm ${OUT}/${target}_treat_pileup.bdg ${OUT}/${target}_control_lambda.bdg ${OUT}/${target}_subtract.bdg

fi

