#!/bin/bash

#get-protocol-strand-interpretation.sh #{{{
usage() { echo "Usage: $0 [-b <bam>] [-g <GTF file>]" 1>&2; exit 1; }

# : means takes an argument but not mandatory (if mandatory will have to check after)
options=':b:g:h'
while getopts $options option
do
    case $option in
        b ) bam=(${OPTARG}) ;;
        g ) gtf=(${OPTARG}) ;;
        h ) usage ;;
        : ) echo "Missing option argument for -$OPTARG" >&2; usage ;;
        \?) echo "Unknown option: -$OPTARG" >&2; usage ;;
        * ) usage ;;
    esac
done

shift $((OPTIND-1))

if [ -z "$bam" ] || [ -z "$gtf" ]
then
    usage
fi

# Check if gene-models file exist
gene_model=${gtf/.gtf/.gene-models.bed12}

if [ ! -f "$gene_model" ]
then
    gtfToGenePred -genePredExt -geneNameAsName2 $gtf ${gene_model/.bed12/.tmp}
    awk '{print $2"\t"$4"\t"$5"\t"$1"\t0\t"$3"\t"$6"\t"$7"\t0\t"$8"\t"$9"\t"$10}' ${gene_model/.bed12/.tmp} >  $gene_model
    rm "${gene_model/.bed12/.tmp}"
fi

infer_experiment.py -r $gene_model -i $bam

