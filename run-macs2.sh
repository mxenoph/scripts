#!/usr/bin/env bash

BROAD=false;
#run-macs2.sh #{{{
usage() { echo "Usage: $0 [-i <IP bam>] [-c <Input bam>] [-p <ncores>] [-g <genome e.g.mm>] [-o <output_path>] [-s <sample name>] [-b <broad-boolean>]" 1>&2; exit 1; }

# : means takes an argument but not mandatory (if mandatory will have to check after)
options=':i:c:g:p:s:o:a:bh'
while getopts $options option
do
    case $option in
        i ) ip=(${OPTARG}) ;;
        c ) CTRL=${OPTARG} ;;
        g ) genome=${OPTARG} ;;
        o ) OUT=${OPTARG} ;;
        b ) BROAD=true ;;
        p ) N_JOBS=${OPTARG} ;;
        s ) SAMPLE=${OPTARG} ;;
        a ) ASSEMBLY=${OPTARG} ;;
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
target="${base%.*}.macs2"
options=''

# If control provided then will pass it to macs#{{{
if [ ! -z "$CTRL" ]
then
    # --bdg will output also bedgraph files for ip and control and --SPMR specifies that the
    # fragment profiles for the pileup are per million reads
    options="$options -c $CTRL --bdg --SPMR"
fi
#}}}

# Controlling output directory#{{{
if [ ! -z "$OUT" ]
then
    if [[ ! $OUT =~ \/macs$ ]]
    then 
        OUT="$(readlink -m ${OUT})/macs"
        mkdir -p $OUT
    fi
    options="$options --outdir $OUT"
else
    options="$options --outdir $PWD"
fi
#}}}

# Controlling peak shape#{{{
if [ ${BROAD} == true ]
then
    options="$options --broad -q 0.05 --broad-cutoff 0.05"
else
    options="$options -q 0.01 --call-summits"
fi
#}}}

#Run macs
macs2 callpeak \
    -t ${ip[*]} \
    --g=$genome -n=$target \
    --verbose 0 $options

if [ ! -z "$CTRL" ] & [ ! -z "$ASSEMBLY" ]
then
    ASSEMBLY=$(echo $ASSEMBLY | tr '[:upper:]' '[:lower:]')
    macs2 bdgcmp -t "${OUT}/${target}_treat_pileup.bdg" -c "${OUT}/${target}_control_lambda.bdg" \
        -m subtract --outdir $OUT --o-prefix ${target}
    
    # Adding track line to broadPeak

    if [ ${BROAD} == true ]
    then
        sed -i "1s/^/track type=broadPeak visibility=3 db=${ASSEMBLY} name=\"${target}_peaks.broadPeak\" description=\"${target} ${ASSEMBLY}\"/" ${OUT}/${target}_peaks.broadPeak
    fi

    chrom_sizes=$(grep ${ASSEMBLY} /nfs/research2/bertone/user/mxenoph/genome_dir/assemblies-annotations.config | cut -f 3)
    bedGraphToBigWig ${OUT}/${target}_subtract.bdg ${chrom_sizes} ${OUT}/${target}_subtract.bw

    rm ${OUT}/${target}_treat_pileup.bdg ${OUT}/${target}_control_lambda.bdg ${OUT}/${target}_subtract.bdg

fi

