#!/usr/bin/env bash

#make-scores-negative.sh #{{{
usage() { echo "Usage: $0 [-b <bigwig>] [-c <chromosomes lengths>]" 1>&2; exit 1; }

# : means takes an argument but not mandatory (if mandatory will have to check after)
options=':b:c:h'
while getopts $options option
do
    case $option in
        b ) bw=(${OPTARG}) ;;
        c ) chr=(${OPTARG}) ;;
        h ) usage ;;
        : ) echo "Missing option argument for -$OPTARG" >&2; usage ;;
        \?) echo "Unknown option: -$OPTARG" >&2; usage ;;
        * ) usage ;;
    esac
done
shift $((OPTIND-1))

# Check that mandatory arguments were provided
if [ -z "$bw" ]
then
    usage
fi

if [ -z "$chr" ]
then
    chr='/nfs/research2/bertone/user/mxenoph/common/genome/MM10/MM10.genome'
fi
#}}}

bigWigToBedGraph ${bw} ${bw/.bw/.bedgraph}
Rscript -e "bw = read.table('${bw/.bw/.bedgraph}', sep = '\t', head=F); bw[,'V4'] = bw[,'V4'] * -1; write.table(bw, file='${bw/.bw/.bedgraph}', quote=FALSE, col.names=FALSE, row.names = FALSE, sep = '\t')"
bedGraphToBigWig  ${bw/.bw/.bedgraph} ${chr} ${bw/.bw/_ng.bw}
rm ${bw/.bw/.bedgraph}

