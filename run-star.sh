#!/usr/bin/env bash

#run-star.sh -f "data/pe/2lox_Epi_*" -g /nfs/research2/bertone/user/mxenoph/common/genome/MM10/ -p 4 -s test -t mm10/#{{{
usage() { echo "Usage: $0 [-f <\"data/pe/2lox_Epi_*\">] [-g <star_genome_dir>] [-p <ncores>] [-s <sample_name>] [-o <output_path>]" 1>&2; exit 1; }

# : means takes an argument but not mandatory (if mandatory will have to check after)
options=':f:g:p:s:o:h'
while getopts $options option
do
    case $option in
        f ) fastq=(${OPTARG}) ;;
        g ) genome=${OPTARG} ;;
        p ) N_JOBS=${OPTARG} ;;
        s ) SAMPLE=${OPTARG} ;;
        o ) TARGET=${OPTARG} ;;
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

options=''

# Functions #{{{

# Checking compressed#{{{
test_compressed(){
   fq="$1"
    
   if [[ $(readlink $fq) =~ ".gz" ]]
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
    if [[ "${fq[0]/_1_sequence*/}" != "${fq[1]/_2_sequence*/}" ]]
    then
        echo 'ERROR: Fastq files provided are either not for mate reads or not in mate_1 mate_2 order.' 1>&2
        exit 1;
    fi
}
#}}}

test_input(){
   fq=("${@}")

    if [[ ${#fq[@]} -eq 2 ]]
    then
        if [ ! -e $(readlink ${fastq[0]}) ] || [ ! -e $(readlink ${fastq[0]}) ]
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
                options="$options --readFilesCommand gzip -c"
                echo "Input files are compressed. Decompressing internally."
            else
                echo "${fastq[1]} uncompressed"
                echo gzip $(readlink ${fastq[1]})
                if [ -f ${fastq[1]} ]
                then
                    echo 'this is symlink, change the symlink after compression to point to the right file'
                    echo ln -s "$(readlink ${fastq[1]}).gz"
                fi
            fi
        else
            if test_compressed "${fastq[1]}"
            then
                echo "${fastq[1]} compressed"
                echo gzip $(readlink ${fastq[0]})
                if [ -f ${fastq[0]} ]
                then
                    echo 'this is symlink, change the symlink after compression to point to the right file'
                    echo ln -s "$(readlink ${fastq[0]}).gz"
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
                options="$options --readFilesCommand gzip -c"
                echo "Input file is compressed. Decompressing internally."
        fi
    fi

}
#}}}


if [ -z "$SAMPLE" ]
then
    SAMPLE='test'
fi

TARGET="${TARGET%%\/}"/star/1pass
#echo ${TARGET}
#mkdir -p "${TARGET}"
#cd "${TARGET}"
output="${TARGET}/${SAMPLE}"

test_input "${fastq[@]}"
#annotation='data/Mus_musculus.GRCm38.75.repeats.gff'

echo STAR --genomeDir $genome --readFilesIn ${fastq[@]} --runThreadN $N_JOBS
