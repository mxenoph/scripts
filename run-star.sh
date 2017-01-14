#!/usr/bin/env bash

#run-star.sh -f "data/pe/2lox_Epi_*" -g /nfs/research2/bertone/user/mxenoph/common/genome/MM10/ -p 4 -s test -o mm10/#{{{
usage() { echo "Usage: $0 [-f <\"data/pe/2lox_Epi_*\">] [-g <star_genome_dir>] [-p <ncores>] [-s <sample_name>] [-o <output_path>] [-m <mode (e.g.1pass, 2pass) default:1pass>]" 1>&2; exit 1; }

# : means takes an argument but not mandatory (if mandatory will have to check after)
options=':f:g:p:s:o:m:h'
while getopts $options option
do
    case $option in
        f ) fastq=(${OPTARG}) ;;
        g ) genome=${OPTARG} ;;
        p ) N_JOBS=${OPTARG} ;;
        s ) SAMPLE=${OPTARG} ;;
        o ) TARGET=${OPTARG} ;;
        m ) MODE=${OPTARG} ;;
        h ) usage ;;
        : ) echo "Missing option argument for -$OPTARG" >&2; usage ;;
        \?) echo "Unknown option: -$OPTARG" >&2; usage ;;
        * ) usage ;;
    esac
done
shift $((OPTIND-1))

# Check that mandatory arguments were provided
if [ -z "$fastq" ] || [ -z "$genome" ]
then
    usage
fi

#}}}

# Functions #{{{

# Checking compressed#{{{
test_compressed(){
   fq="$1"
    
   if [[ $(realpath $fq) =~ ".gz" ]]
   then
       # 0 values is true in bash
       return 0
   else
       return 1
   fi
}
#}}}

# Checking pe format#{{{
test_pe(){
   fq=("${@}")
    
    # Sanity check that the fastq provided are for mate reads
    if [[ "${fq[0]/_1*/}" != "${fq[1]/_2*/}" ]]
    then
        echo 'ERROR: Fastq files provided are either not for mate reads or not in mate_1 mate_2 order.' 1>&2
        exit 1;
    fi
}
#}}}

# testing paired end validity and compression#{{{
test_input(){
   fq=("${@}")

    if [[ ${#fq[@]} -eq 2 ]]
    then
        if [ ! -e $(realpath ${fastq[0]}) ] || [ ! -e $(realpath ${fastq[1]}) ]
        then
            echo "ERROR: One or both the provided fastq files do not exist. Please check the path and the validity of symbolic links if used." 1>&2
            exit;
        fi

        test_pe "${fq[@]}"

        if test_compressed "${fastq[0]}"
        then
            if test_compressed "${fastq[1]}"
            then
                # if file compressed (only checking for gz in the function) unzip keeping original files and sending to stdout
                options="$options --readFilesCommand zcat"
                echo "Input files are compressed. Decompressing internally."
            else
                echo "$(realpath ${fastq[1]}) is uncompressed (${fastq[0]} is compressed). Compressing now."
                gzip $(realpath ${fastq[1]})
                
                # check if file is symbolic link with -h
                if [ -h ${fastq[1]} ]
                then
                    echo -e 'After compression change the symlink to point to the compressed file.'
                    # updating the already existing symbolic link so that it points to the gzip file
                    ln -sfn "$(realpath ${fastq[1]}).gz" ${fastq[1]}
                    options="$options --readFilesCommand zcat"
                fi
            fi
        else
            if test_compressed "${fastq[1]}"
            then
                echo "$(realpath ${fastq[0]}) is uncompressed (${fastq[1]} is compressed). Compressing now."
                gzip $(realpath ${fastq[0]})

                if [ -h ${fastq[0]} ]
                then
                    echo -e 'After compression change the symlink to point to the compressed file.'
                    ln -sfn "$(realpath ${fastq[0]}).gz ${fastq[0]}"
                    options="$options --readFilesCommand zcat"
                fi
            fi
        fi
    elif [[ ${#fq[@]} -gt 2 ]]
    then
        echo "ERROR: $0 can only run on fastq file at a time." 1>&2
        exit
    else
        if test_compressed ${fastq[0]}
        then
                options="$options --readFilesCommand zcat"
                echo "Input file is compressed. Decompressing internally."
        fi
    fi

}
#}}}

# First pass with annotation because the index was generated with annotation #{{{
first_pass(){
    fq=("${@}")
    
    TARGET="${TARGET%%\/}"/star/1pass
    mkdir -p ${TARGET}
    output="${TARGET}/${SAMPLE}"

    # encode uses --outFilterType BySJout but that will keep only those reads that contain junctions that passed filtering into SJ.out.tab
    # does that mean it will throw away reads not containing a junction?
    encode_options="--outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --outFileNamePrefix ${output}"

    if [[ ${#fq[@]} -eq 2 ]]
    then
        encode_options="${encode_options} --alignMatesGapMax 1000000 --outSAMtype BAM Unsorted"
    else
        encode_options="${encode_options} --outSAMtype BAM SortedByCoordinate"
    fi

    STAR --genomeDir $genome --readFilesIn $(realpath ${fq[@]}) --runThreadN $N_JOBS ${options} ${encode_options}
}
#}}}

# Per-sample 2pass #{{{
second_pass(){
    fq=("${@}")
    
    TARGET="${TARGET%%\/}"/star/2pass
    mkdir -p ${TARGET}
    output="${TARGET}/${SAMPLE}"

    # encode uses --outFilterType BySJout but that will keep only those reads that contain junctions that passed filtering into SJ.out.tab
    # does that mean it will throw away reads not containing a junction?
    encode_options="--outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --outFileNamePrefix ${output} --twopassMode Basic"

    if [[ ${#fq[@]} -eq 2 ]]
    then
        encode_options="${encode_options} --alignMatesGapMax 1000000 --outSAMtype BAM Unsorted"
    else
        encode_options="${encode_options} --outSAMtype BAM SortedByCoordinate"
    fi

    STAR --genomeDir $genome --readFilesIn $(realpath ${fq[@]}) --runThreadN $N_JOBS ${options} ${encode_options}
}
#}}}

#}}}


if [ -z "$SAMPLE" ]
then
    SAMPLE=$(basename ${fastq[0]})
    SAMPLE=${SAMPLE/.*/}
fi

if [ -z "$MODE" ]
then
    MODE='2pass'
fi

if [ -z "$TARGET" ]
then
    TARGET=$PWD
fi

# It will set the decompression internally if the files are compressed
options=''
test_input "${fastq[@]}"

if [[ "$MODE" == "1pass" ]]
then
    echo first pass mode
    first_pass "${fastq[@]}"
elif [[ "$MODE" == "2pass" ]]
then
    echo second pass mode
    second_pass "${fastq[@]}"
else
    echo "Unknown mode provided."
fi

