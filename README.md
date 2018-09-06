bctools
=======

A toolbox for correcting barcodes in 10X linked-read sequencing data.




Prerequisites
-------------

* GCC version >= 4.9 (supports C++14)
* SeqAn core library, version 2.3.1? (https://github.com/seqan/seqan)
* SDSL - Succinct Data Structure Lirbrary (https://github.com/simongog/sdsl-lite)
* kseq.h from HTSlib (https://github.com/samtools/htslib)




Installation
------------

1. Download the Seqan core library. You do not need to follow the SeqAn install instructions. You only need the directory .../include/seqan with all its content (the SeqAn core library).
2. Download and install the SDSL.
3. Download HTSlib or just put the kseq.h header file into a folder named htslib.
3. Edit lines 14-17 in the Makefile to point to the directories of SeqAn, SDSL and HTSlib.
4. Run 'make' in the bctools directory.

If everything is setup correctly, this will create the binary 'bctools'.




Usage
-----

The only input needed for barcode correction is a pair of barcoded FASTQ files generated on the 10X Chromium platform.
Optionally, you can specify a barcode whitelist file.

The program consists of several commands, which are listed when running

    ./bctools --help

For a short description of each command and an overview of arguments and options, you can run

    ./bctools <COMMAND> --help

If you need the output to be sorted and/or converted to SAM, BAM, or (gzipped) FASTQ format, you can run the provided bash script. For a short description of options and arguments of this script run

    ./scripts/run_bctools -h

### The whitelist command

    ./bctools whitelist [OPTIONS] <FASTQ 1 file>

Creates a barcode whitelist based on barcode occurence in the data.
Creating a whitelist from your data is recommended (rather than using the 10X whitelist) to reduce the number of alternatives during correction and prevents false corrections.

### The index command

    ./bctools index [OPTIONS] <whitelist file>

Creates a barcode index from the given barcode whitelist and writes it to disk. This command is optional as the index can be created on the fly in the 'correct' command.

### The correct command

    ./bctools correct [OPTIONS] <whitelist file> <FASTQ 1 file> <FASTQ 2 file>

Corrects barcodes of the given barcoded read pair data using the specified barcode whitelist. A barcode index is computed on the fly unless index files are present for the specified barcode whitelist. The output is a tab-separated file holding one read pair per line as decribed below.

### The stats command

    ./bctools stats [OPTIONS] <Corrected (gzipped) FASTQ 1 file>
    ./bctools stats [OPTIONS] <Corrected SAM/BAM file>
    ./bctools stats [OPTIONS] <Corrected TSV file>

Computes the number of read pairs with whitelisted, corrected and unrecognized barcodes, a barcode occurrence histogram and counts quality values of corrected barcode positions.




Example
-------

    mkdir bctools_example && cd bctools_example/
    ln -s /path/to/first.fq.gz
    ln -s /path/to/second.fq.gz

    ./bctools whitelist -o whitelist.txt first.fq.gz
    ./bctools correct whitelist.txt first.fq.gz second.fq.gz > corrected.tsv

Using the bash script to create a BAM file sorted by the corrected barcode sequence:

    ./script/run_bctools -f bam first.fq.gz second.fq.gz


Output format
-------------

The output format of the correct command is a simple tab-separated format, where each read pair and its barcode information is given on a single line.
The fields are as follows:

Field | Description
--- | ---
READ NAME | The read or query name taken from the FASTQ file and cropped at the first whitespace. 
CORRECTED BARCODE | A comma separated list of possible barcode corrections. If the raw barcode is whitelisted, the value of this field is identical to the RAW BARCODE field. An asterisk '*' indicates that the barcode is not whitelisted and correction was unsuccessful.
RAW BARCODE | The first 16 base pairs of the first read in the read pair.
7-MER SPACER | The seven base pairs following the first 16 base pairs of the first read in the read pair.
TRIMMED FIRST READ | The remaining base pairs of the first read in the read pair after trimming the barcode and 7-mer spacer sequence.
SECOND READ | The second read sequence.
BARCODE QUALITY STRING | The first 16 values of the quality string of the first read in the read pair.
7-MER SPACER QUALITY STRING |  The seven values following the first 16 values of the quality string of the first read in the read pair.
TRIMMED FIRST READ QUALITY STRING | The remaining quality string after trimming the barcode and 7-mer spacer quality values.
SECOND READ QUALITY STRING | The quality string of the second read in the read pair.


Contact
-------

For questions and comments contact birte.kehr [at] bihealth.de or create an issue.
