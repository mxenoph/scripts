#! /usr/bin/env bash

# For creating pseudoreplicates from file

#make-pseudoreps.sh #{{{
usage() { echo "Usage: $0 [-b <bam file to be split>] [-o <output directory>] [-n <number of pseudoreplicates>]" 1>&2; exit 1; }

# : means takes an argument but not mandatory (if mandatory will have to check after)
options=':b:o:n:h'
while getopts $options option
do
    case $option in
        b ) bam=(${OPTARG}) ;;
        o ) output_path=${OPTARG} ;;
        n ) NUM=${OPTARG} ;;
        h ) usage ;;
        : ) echo "Missing option argument for -$OPTARG" >&2; usage ;;
        \?) echo "Unknown option: -$OPTARG" >&2; usage ;;
        * ) usage ;;
    esac
done
shift $((OPTIND-1))

if [ -z "$bam" ] || [ -z "$output_path" ]
then
    usage
fi

if [ -z "$NUM" ]
then
    NUM=2
fi

#}}}

target=$(basename ${bam/.bam/})

# This will shuffle the lines in the file and split it into two parts
total_reads=$(samtools view -c ${bam})
subset_reads=$(( total_reads / ${NUM} ))
command="samtools view ${bam} | shuf | split -d -l ${subset_reads} - ${output_path}/${target}"
eval ${command}
samtools view -H ${bam} > ${output_path}/${target}.sam_header
# In the last sed command using @ as the delimiter instead of / because $command (which is expanded in the substitution) contains / and so sed thinks the command ends before it actually does
tail -n 1 ${output_path}/${target}.sam_header | grep '@PG' | sed 's/ID:.*\tVN/ID:shuf-split\tVN/g' | sed 's/VN:.*CL/VN:unk\tCL/' | sed "s@CL:.*@CL:\"$command\"@g" > ${output_path}/${target}.sam_PG_line

cat ${output_path}/${target}.sam_header ${output_path}/${target}.sam_PG_line ${output_path}/${target}00 > ${output_path}/${target}00.sam
cat ${output_path}/${target}.sam_header ${output_path}/${target}.sam_PG_line ${output_path}/${target}01 > ${output_path}/${target}01.sam

samtools view -bSo ${output_path}/${target}00.bam ${output_path}/${target}00.sam
samtools view -bSo ${output_path}/${target}01.bam ${output_path}/${target}01.sam

#rm ${output_path}/${target}00.sam
#rm ${output_path}/${target}01.sam
#rm ${output_path}/${target}.sam_PG_line
#rm ${output_path}/${target}.sam_header

