#!/usr/bin/env bash

#run-bowtie.sh #{{{
usage() { echo "Usage: $0 [-b <bowtie folder>]" 1>&2; exit 1; }

# : means takes an argument but not mandatory (if mandatory will have to check after)
options=':b:h'
while getopts $options option
do
    case $option in
        b ) bowtie=(${OPTARG}) ;;
        h ) usage ;;
        : ) echo "Missing option argument for -$OPTARG" >&2; usage ;;
        \?) echo "Unknown option: -$OPTARG" >&2; usage ;;
        * ) usage ;;
    esac
done
shift $((OPTIND-1))

if [ -z "$bowtie" ]
then
    usage
fi

if [ ! -d "$bowtie" ]
then
    usage
fi

# }}}

#get the name of the script for the logs
output=${bowtie}/"bowtie-alignment-statistics.tsv"
echo -e "filename\tall\tmapped\tfail\tmore_mult\tdup_removed" > ${output}
log_files=$(ls ${bowtie}/*.stderr)
echo ${log_files}

for file in ${log_files[@]}
do
	filename=`basename $file`

	# For bowtie result
    index=0
	for stat in `awk '$0 ~ /^# reads/ { print $0 }' $file | cut -d ':' -f 2 | cut -d ' ' -f 2`;
	do
		per_file_stats[${index}]=${stat}
		index=$((index+1))
	done

	picard_file=${bowtie}/ddup/${filename/.stderr/_ddup.bam}
	echo $picard_file
    # This counts the number of alignments in the files. Because the picard file was generated by
    # run-bowtie.sh which removes duplicates, the count = mapped duplicate free reads + unmapped reads
    number_ddup_reads=$(samtools view -c ${picard_file})

    # per_file_stats[0:total reads, 1: unique, 2:failed, 3:suppressed]
	echo -e "${filename}\t${per_file_stats[0]}\t${per_file_stats[1]}\t${per_file_stats[2]}\t${per_file_stats[3]}\t$number_ddup_reads" >> $output
done;
