#!/bin/bash

set -e
set -u
set -o pipefail

################################################################################

usage()
{
    echo "Barcode correction, sorting, file conversion" 1>&2
    echo "============================================" 1>&2
    echo "" 1>&2
    echo "SYNOPSIS" 1>&2
    echo "    ./$0 [OPTIONS] <FASTQ1> <FASTQ2>" 1>&2
    echo "" 1>&2
    echo "    Script for correcting barcodes with bctools, sorting the output by barcode" 1>&2
    echo "    with Unix sort, and file conversion to (gzipped) FASTQ, SAM, or BAM." 1>&2
    echo "" 1>&2
    echo "    -h" 1>&2
    echo "        Display this help message." 1>&2
    echo "" 1>&2
    echo "BARCODE CORRECTION OPTIONS" 1>&2
    echo "    -b  STR" 1>&2
    echo "        Path to bctools program, e.g. '/path/to/bctools'. Default: bctools" 1>&2
    echo "    -w  STR" 1>&2
    echo "        Name of whitelist file. If file does not exist, whitelist is inferred" 1>&2
    echo "        from the data and written to this file." 1>&2
    echo "        Default: <input prefix>.whitelist.txt" 1>&2
    echo "    -a  NUM" 1>&2
    echo "        Maximum number of alternative barcode corrections." 1>&2
    echo "        In range [1..16]. Default: 4." 1>&2
    echo "" 1>&2
    echo "SORTING OPTIONS" 1>&2
    echo "    -n" 1>&2
    echo "        Do not sort output of bctools by barcode." 1>&2
    echo "    -S  SIZE" 1>&2
    echo "        Size for main memory buffer, e.g. '8G'. Value is passed to Unix sort" 1>&2
    echo "        -S, --buffer-size option. Default: 4G" 1>&2
    echo "    -T  DIR" 1>&2
    echo "        Directory for tempoarary files, e.g. '/path/to/tmp/'. Value is passed to" 1>&2
    echo "        Unix sort -T, --temporary-directory option." 1>&2
    echo "" 1>&2
    echo "OUTPUT OPTIONS" 1>&2
    echo "    -o  STR" 1>&2
    echo "        Output prefix. Default: <input prefix>.corrected" 1>&2
    echo "    -f  STR" 1>&2
    echo "        Ouptut format. One of: fastq, fastq.gz, sam, bam, tsv" 1>&2
    echo "    -s  STR" 1>&2
    echo "        Path to samtools program for BAM conversion, e.g. '/path/to/samtools'." 1>&2
    echo "        Necessary only if output format is 'bam'. Default: samtools" 1>&2
    echo "" 1>&2
    echo "INFO" 1>&2
    echo "    Created on Aug 31, 2018" 1>&2
    exit 1
}

################################################################################
# FASTQ conversion with compression

convert_and_compress_first()
{
    awk '{print "@" $1 " BX:Z:" $2 " RX:Z:" $3 " QX:Z:" $7 " TR:Z:" $4 " TQ:Z:" $8 "\n" $5 "\n+\n" $9}' \
    | gzip -c \
    > ${outprefix}.1.fastq.gz
}
export -f convert_and_compress_first

convert_and_compress_second()
{
    awk '{print "@" $1 "\n" $6 "\n+\n" $10}' \
    | gzip -c \
    > ${outprefix}.2.fastq.gz
}
export -f convert_and_compress_second

################################################################################
# FASTQ conversion without compression

convert_first()
{
    awk '{print "@" $1 " BX:Z:" $2 " RX:Z:" $3 " QX:Z:" $7 " TR:Z:" $4 " TQ:Z:" $8 "\n" $5 "\n+\n" $9}' \
    > ${outprefix}.1.fastq
}
export -f convert_first

convert_second()
{
    awk '{print "@" $1 "\n" $6 "\n+\n" $10}' \
    > ${outprefix}.2.fastq
}
export -f convert_second

################################################################################
# SAM conversion

convert_sam()
{
    awk -v prefix=${outprefix} '
       BEGIN {OFS="\t";
           print "@HD", "VN:1.3", "SO:unknown";
           print "@RG", "ID:" name ":" prefix, "SM:" prefix, "LB:" prefix, "PU:1", "PL:ILLUMINA";
       } {
           bx="";
           gsub(",", "-1,", $2)
           if ($2 != "*") bx="BX:Z:" $2 "-1\t";
           print $1, "68",  "*", "0", "0", length($5) "S", "*", "0", "-1", $5, $9,  "RG:Z:" prefix, "TR:Z:" $4, "TQ:Z:" $8, bx "RX:Z:" $3, "QX:Z:" $7;
           print $1, "132", "*", "0", "0", length($6) "S", "*", "0", "-1", $6, $10, "RG:Z:" prefix,                         bx "RX:Z:" $3, "QX:Z:" $7;
       }'
}

################################################################################
### MAIN ###

bctools="bctools"
whitelist="-"
cutoff="inferred"
alts=4
sort="on"
buffersize=4G
tempdir="-"
outprefix="-"
format="fastq.gz"
samtools="samtools"

while getopts 'hb:w:c:a:nS:T:o:f:s:' OPTION; do
    case "${OPTION}" in
        b)
            bctools="${OPTARG}"
            ;;
        w)
            whitelist="${OPTARG}"
            ;;
        c)
            cutoff="${OPTARG}"
            ;;
        a)
            alts="${OPTARG}"
            ;;
        n)
            sort="off"
            ;;
        S)
            buffersize="${OPTARG}"
            ;;
        T)
            tempdir="${OPTARG}"
            ;;
        o)
            outprefix="${OPTARG}"
            ;;
        f)
            format="${OPTARG}"
            ;;
        s)
            samtools="${OPTARG}"
            ;;
        *)
            usage
            ;;
    esac
done
shift "$((${OPTIND} -1))"
if [ $# -ne 2 ]; then
    echo "ERROR: Please specify two input FASTQ files (options go first)" 1>&2
    echo "       or -h to display a help message." 1>&2
    exit 1
fi

FASTQ1=$1
FASTQ2=$2

# Set outprefix to '<input prefix>corrected'.
if [ "${outprefix}" = "-" ]; then
    outprefix=`printf "%s\n%s\n" "${FASTQ1}" "${FASTQ2}" | sed -e 'N;s/^\(.*\).*\n\1.*$/\1/'`
    outprefix="${outprefix}corrected"
fi

# Set whitelist filename to '<input prefix>whitelist.txt'.
if [ "${whitelist}" = "-" ]; then
    whitelist=`printf "%s\n%s\n" "${FASTQ1}" "${FASTQ2}" | sed -e 'N;s/^\(.*\).*\n\1.*$/\1/'`
    whitelist="${whitelist}whitelist.txt"
fi

# Set the sorting options.
sortopt=""
if [ "${buffersize}" != '-' ]; then
    sortopt="${sortopt} --buffer-size=${buffersize}"
fi
if [ "${tempdir}" != '-' ]; then
    sortopt="${sortopt} --temporary-directory=${tempdir}"
fi

echo "Barcode correction, sorting, file conversion" 1>&2
echo "============================================" 1>&2
echo "" 1>&2
echo "Arguments:" 1>&2
echo "    FASTQ1: $FASTQ1" 1>&2
echo "    FASTQ2: $FASTQ2" 1>&2
echo "" 1>&2
echo "Options:" 1>&2
echo "    Path to bctools program:       $bctools" 1>&2
echo "    Whitelist file:                $whitelist" 1>&2
echo "    Whitelist cutoff:              $cutoff" 1>&2
echo "    Max. num. of alt. corrections: $alts" 1>&2
echo "    Sorting:                       $sort" 1>&2
echo "    Sorting buffer size:           $buffersize" 1>&2
echo "    Sorting tmp directory:         $tempdir" 1>&2
echo "    Output prefix:                 $outprefix" 1>&2
echo "    Output format:                 $format" 1>&2
echo "    Path to samtools program:      $samtools" 1>&2
echo "" 1>&2

# Create a barcode whitelist by running bctools whitelist.
if [ ! -f "${whitelist}" ]; then
    if [ "${cutoff}" = "inferred" ]; then
        ${bctools} whitelist -o ${whitelist} ${FASTQ1}
    else
        ${bctools} whitelist -o ${whitelist} -c ${cutoff} ${FASTQ1}
    fi
else
    echo "NOTE: Whitelist file exists. Skipping '${bctools} whitelist' command." 1>&2
    echo "      If you want the whitelist to be created, remove the existing file" 1>&2
    echo "      or specify another filename." 1>&2
fi
echo "" 1>&2

# Run barcode correction, sorting and file conversion as specified.
if [ ${sort} = "on" ]; then
    case "${format}" in
        "fastq.gz")
            ${bctools} correct -a ${alts} ${whitelist} ${FASTQ1} ${FASTQ2} \
            | sort -k3,3 ${sortopt} \
            | tee >(convert_and_compress_first) \
            | convert_and_compress_second
            echo "" 1>&2
            echo "Output written to '${outprefix}.1.fastq.gz' and '${outprefix}.2.fastq.gz'." 1>&2
            ;;
        "fastq")
            ${bctools} correct -a ${alts} ${whitelist} ${FASTQ1} ${FASTQ2} \
            | sort -k3,3 ${sortopt} \
            | tee >(convert_first) \
            | convert_second
            echo "" 1>&2
            echo "Output written to '${outprefix}.1.fastq' and '${outprefix}.2.fastq'." 1>&2
            ;;
        "sam")
            ${bctools} correct -a ${alts} ${whitelist} ${FASTQ1} ${FASTQ2} \
            | sort -k3,3 ${sortopt} \
            | convert_sam \
            > ${outprefix}.sam
            echo "" 1>&2
            echo "Output written to '${outprefix}.sam'." 1>&2
            ;;
        "bam")
            ${bctools} correct -a ${alts} ${whitelist} ${FASTQ1} ${FASTQ2} \
            | sort -k3,3 ${sortopt} \
            | convert_sam \
            | ${samtools} view -Sb \
            > ${outprefix}.bam
            echo "" 1>&2
            echo "Output written to '${outprefix}.bam'." 1>&2
            ;;
        "tsv")
            ${bctools} correct -a ${alts} ${whitelist} ${FASTQ1} ${FASTQ2} \
            | sort -k3,3 ${sortopt} \
            > ${outprefix}.tsv
            echo "" 1>&2
            echo "Output written to '${outprefix}.tsv'." 1>&2
    esac
else
    case "${format}" in
        "fastq.gz")
            ${bctools} correct -a ${alts} ${whitelist} ${FASTQ1} ${FASTQ2} \
            | tee >(convert_and_compress_first) \
            | convert_and_compress_second
            echo "" 1>&2
            echo "Output written to '${outprefix}.1.fastq.gz' and '${outprefix}.2.fastq.gz'." 1>&2
            ;;
        "fastq")
            ${bctools} correct -a ${alts} ${whitelist} ${FASTQ1} ${FASTQ2} \
            | tee >(convert_first) \
            | convert_second
            echo "" 1>&2
            echo "Output written to '${outprefix}.1.fastq' and '${outprefix}.2.fastq'." 1>&2
            ;;
        "sam")
            ${bctools} correct -a ${alts} ${whitelist} ${FASTQ1} ${FASTQ2} \
            | convert_sam \
            > ${outprefix}.sam
            echo "" 1>&2
            echo "Output written to '${outprefix}.sam'." 1>&2
            ;;
        "bam")
            ${bctools} correct -a ${alts} ${whitelist} ${FASTQ1} ${FASTQ2} \
            | convert_sam \
            | samtools view -Sb \
            > ${outprefix}.bam
            echo "" 1>&2
            echo "Output written to '${outprefix}.bam'." 1>&2
            ;;
        "tsv")
            ${bctools} correct -a ${alts} ${whitelist} ${FASTQ1} ${FASTQ2} \
            > ${outprefix}.tsv
            echo "" 1>&2
            echo "Output written to '${outprefix}.tsv'." 1>&2
    esac
fi

