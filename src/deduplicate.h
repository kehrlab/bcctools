#ifndef DEDUPLICATE_H_
#define DEDUPLICATE_H_

#include <cassert>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/align.h>
#include "utils.h"

struct DupPrefix {};
struct DupReadName {};
struct DupBoth {};

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

bool isCandidateDup(ReadPair & a, ReadPair & b, unsigned minMatches, unsigned maxOffset, DupPrefix const);
bool isCandidateDup(ReadPair & a, ReadPair & b, unsigned minMatches, unsigned maxOffset, DupReadName const);
bool isCandidateDup(ReadPair & a, ReadPair & b, unsigned minMatches, unsigned maxOffset, DupBoth const);
int getQual(ReadPair & rp, unsigned minQual);

template <typename TDupTag>
void find_duplicates(std::vector<ReadPair> & readPairs, unsigned minMatches, unsigned maxOffset, double maxDiffRate, unsigned minQual, TDupTag const tag)
{
    std::sort(std::begin(readPairs), std::end(readPairs), compareBySeq);

    // Align pairs of read pairs if the first 'minMatches' bases of the first reads in the two pairs match.
    for (unsigned i = 0; i < readPairs.size(); ++i)
    {
        if (readPairs[i].isDup)
            continue;

        for (unsigned j = i+1; j < readPairs.size() && isCandidateDup(readPairs[i], readPairs[j], minMatches, maxOffset, tag); ++j)
        {
            if (readPairs[j].isDup)
                continue;

            // Align the two read pairs
            int score1 = seqan::globalAlignmentScore(readPairs[i].read1, readPairs[j].read1, seqan::Score<int, seqan::Simple>(0, -1, -1), -3, 3);

            // // Debug output.
            // int s2 = seqan::globalAlignmentScore(readPairs[i].read2, readPairs[j].read2, seqan::Score<int, seqan::Simple>(0, -1, -1), -3, 3);
            // if (score1 > -50 && s2 > -50)
            // {
            //     std::cout << readPairs[i].qname << " " << readPairs[i].read1 << " " << readPairs[i].read2 << std::endl;
            //     std::cout << readPairs[j].qname << " " << readPairs[j].read1 << " " << readPairs[j].read2 << std::endl << std::endl;
            //     std::cout << readPairs[i].qname << " " << readPairs[i].qual1 << " " << readPairs[i].qual2 << std::endl;
            //     std::cout << readPairs[j].qname << " " << readPairs[j].qual1 << " " << readPairs[j].qual2 << std::endl;
            //     std::cout << score1 << "  " << s2 << std::endl << std::endl;
            // }

            if (score1 > -maxDiffRate * length(readPairs[i].read1))
            {
                int score2 = seqan::globalAlignmentScore(readPairs[i].read2, readPairs[j].read2, seqan::Score<int, seqan::Simple>(0, -1, -1), -3, 3);
                if (score2 > -maxDiffRate * length(readPairs[i].read2))
                {
                    // Set lower quality read pair as duplicate
                    if (getQual(readPairs[i], minQual) < getQual(readPairs[j], minQual))
                    {
                        readPairs[i].isDup = true;
                        break; // the loop over j with constant i
                    }
                    else
                        readPairs[j].isDup = true;
                }
            }
        }
    }
};

#endif // DEDUPLICATE_H_;