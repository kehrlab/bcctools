#include <seqan/align.h>

#include "deduplicate.h"

bool compareBySeq(ReadPair const& first, ReadPair const& second)
{
    if (first.read1 != second.read1)
        return first.read1 < second.read1;
    return first.read2 < second.read2;
}

bool compareByName(ReadPair const& first, ReadPair const& second)
{
    if (first.qnameSplit.size() != second.qnameSplit.size())
        return first.qname < second.qname;

    for (size_t i = 0; i < first.qnameSplit.size(); ++i)
    {
        if (first.qnameSplit[i] != second.qnameSplit[i])
        {
            try
            {
                return std::stoi(toCString(first.qnameSplit[i]), nullptr) < std::stoi(toCString(second.qnameSplit[i]), nullptr);
            }
            catch (std::invalid_argument const &)
            {
                return first.qnameSplit[i] < first.qnameSplit[i];
            }
        }
    }
    return false;
}

bool goNext(Tsv_iterator & tsv_it)
{
    if (tsv_it.atEnd)
        return false;

    if (tsv_it.nextPair.cBarcode.size() == 0)
        tsv_it.barcode = "";
    else
        tsv_it.barcode = tsv_it.nextPair.cBarcode[0];

    tsv_it.readPairs.clear();
    tsv_it.readPairs.push_back(tsv_it.nextPair);

    ReadPair nextPair;
    while(readTsvLine(nextPair, tsv_it.in))
    {
        if (tsv_it.barcode == "" || tsv_it.barcode != nextPair.cBarcode[0])
        {
            tsv_it.nextPair = nextPair;
            return true;
        }
        tsv_it.readPairs.push_back(nextPair);
    }

    tsv_it.atEnd = true;
    return true;
}

bool readTsvLine(ReadPair & nextPair, std::ifstream & in)
{
    std::string qname, cBarcode, rBarcode, spacer, read1, read2, qBarcode, qSpacer, qual1, qual2;
    if (in >> qname >> cBarcode >> rBarcode >> spacer >> read1 >> read2 >> qBarcode >> qSpacer >> qual1 >> qual2)
    {
        nextPair = ReadPair(qname, cBarcode, rBarcode, spacer, read1, read2, qBarcode, qSpacer, qual1, qual2);
        return true;
    }
    else
    {
        return false;
    }
}

int getQual(ReadPair & rp, unsigned minQual)
{
    if (rp.avgQual == -1)
    {
        rp.avgQual = 0;

        // First read in pair.
        for (unsigned i = 0; i < length(rp.qual1); ++i)
        {
            // Add qual value if larger than minQual.
            if (rp.qual1[i] > minQual + 33)
                rp.avgQual += rp.qual1[i] - 33;
        }

        // Second read in pair.
        for (unsigned i = 0; i < length(rp.qual2); ++i)
        {
            // Add qual value if larger than minQual.
            if (rp.qual2[i] > minQual + 33)
                rp.avgQual += rp.qual2[i] - 33;
        }
    }

    return rp.avgQual;
}

void find_optical_duplicates(std::vector<ReadPair> & readPairs)
{
    std::sort(std::begin(readPairs), std::end(readPairs), compareByName);
    // TODO
}

void find_sequence_duplicates(std::vector<ReadPair> & readPairs, unsigned minMatches, double maxDiffRate, unsigned minQual)
{
    std::sort(std::begin(readPairs), std::end(readPairs), compareBySeq);

    // Align pairs of read pairs if the first 'minMatches' bases of the first reads in the two pairs match.
    for (unsigned i = 0; i < readPairs.size(); ++i)
    {
        if (readPairs[i].isDup)
            continue;

        for (unsigned j = i+1; j < readPairs.size() && prefix(readPairs[i].read1, minMatches) == prefix(readPairs[j].read1, minMatches); ++j)
        {
            if (readPairs[j].isDup)
                continue;

            // Align the two read pairs
            int score1 = seqan::globalAlignmentScore(readPairs[i].read1, readPairs[j].read1, seqan::Score<int, seqan::Simple>(0, -1, -1), -3, 3);

            // // Debug output.
            // int s2 = seqan::globalAlignmentScore(readPairs[i].read2, readPairs[j].read2, seqan::Score<int, seqan::Simple>(0, -1, -1), -3, 3);
            // if (score1 > -50 && s2 > -50)
            // {
            //     std::cout << readPairs[j].read1 << " " << readPairs[j].read2 << std::endl;
            //     std::cout << readPairs[i].read1 << " " << readPairs[i].read2 << std::endl << std::endl;
            //     std::cout << readPairs[j].qual1 << " " << readPairs[j].qual2 << std::endl;
            //     std::cout << readPairs[i].qual1 << " " << readPairs[i].qual2 << std::endl;
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
}