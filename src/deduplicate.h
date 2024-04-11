#ifndef DEDUPLICATE_H_
#define DEDUPLICATE_H_

#include <cassert>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include "utils.h"

struct ReadPair {
    seqan::CharString qname;
    std::vector<seqan::CharString> qnameSplit;
    std::vector<seqan::DnaString> cBarcode;
    seqan::Dna5String rBarcode;
    seqan::Dna5String spacer;
    seqan::Dna5String read1;
    seqan::Dna5String read2;
    seqan::CharString qBarcode;
    seqan::CharString qSpacer;
    seqan::CharString qual1;
    seqan::CharString qual2;
    bool isDup;
    int avgQual;

    ReadPair () {}

    ReadPair(std::string & qname, std::string & cBarcode, std::string & rBarcode, std::string & spacer, std::string & read1, std::string & read2, std::string & qBarcode, std::string & qSpacer, std::string & qual1, std::string & qual2) :
        qname(qname), rBarcode(rBarcode), spacer(spacer), read1(read1), read2(read2), qBarcode(qBarcode), qSpacer(qSpacer), qual1(qual1), qual2(qual2), isDup(false), avgQual(-1)
    {
        if (seqan::length(cBarcode) > 0)
            this->cBarcode = parseBarcodeList(seqan::toCString(cBarcode));

        seqan::strSplit(qnameSplit, qname, seqan::EqualsChar<':'>());
    }
};


bool compareBySeq(ReadPair const& first, ReadPair const& second);
bool compareByName(ReadPair const& first, ReadPair const& second);

bool readTsvLine(ReadPair & nextPair, std::ifstream & in);

struct Tsv_iterator {
    std::ifstream in;
    ReadPair nextPair;
    bool atEnd;

    seqan::DnaString barcode;
    std::vector<ReadPair> readPairs;
    
    Tsv_iterator(seqan::CharString & filename) :
        in(std::ifstream(toCString(filename)))
    {
        atEnd = !readTsvLine(nextPair, in);
    }
};

bool goNext(Tsv_iterator & tsv_it);

void find_optical_duplicates(std::vector<ReadPair> & readPairs);
void find_sequence_duplicates(std::vector<ReadPair> & readPairs, unsigned minMatches, double maxDiffRate, unsigned minQual);

#endif // DEDUPLICATE_H_;