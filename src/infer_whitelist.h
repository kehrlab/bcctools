#ifndef INFER_WHITELIST_H_
#define INFER_WHITELIST_H_

#include <seqan/sequence.h>
#include <zlib.h>
#include <htslib/kseq.h>

#ifndef KSEQ_GZ
#define KSEQ_GZ
KSEQ_INIT(gzFile, gzread)
#endif  // KSEQ_GZ

void make_histograms(std::vector<unsigned> & allHist,std::vector<unsigned> & wlHist, std::vector<uint16_t> & counts_per_barcode, seqan::CharString & whitelistFile, unsigned bcLength);
double entropy(seqan::DnaString & bc);
unsigned infer_cutoff(std::vector<unsigned> & allHist, std::vector<unsigned> & wlHist);
#endif  // INFER_WHITELIST_H_
