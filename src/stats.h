#ifndef STATS_H_
#define STATS_H_

#include <map>
#include <seqan/sequence.h>
#include <zlib.h>
#include <htslib/kseq.h>

#ifndef KSEQ_GZ
#define KSEQ_GZ
KSEQ_INIT(gzFile, gzread)
#endif  // KSEQ_GZ

struct BarcodeStats
{
    unsigned error_free;
    unsigned one_error;
    unsigned unrecognized;

    seqan::DnaString prev_corrected;
    unsigned prev_count;

    std::vector<unsigned> count_hist;
    std::map<char, std::vector<unsigned> > one_error_hist;

    BarcodeStats() :
        error_free(0), one_error(0), unrecognized(0), prev_count(0)
    {}
};

void write_stats(seqan::CharString & outputFile, BarcodeStats & stats);
void stats_tsv(BarcodeStats & stats, seqan::CharString & inputFile);
void stats_fastq(BarcodeStats & stats, seqan::CharString & inputFile);
void stats_bam(BarcodeStats & stats, seqan::CharString & inputFile);

#endif  // STATS_H_
