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
            if ((unsigned)rp.qual1[i] > minQual + 33)
                rp.avgQual += rp.qual1[i] - 33;
        }

        // Second read in pair.
        for (unsigned i = 0; i < length(rp.qual2); ++i)
        {
            // Add qual value if larger than minQual.
            if ((unsigned)rp.qual2[i] > minQual + 33)
                rp.avgQual += rp.qual2[i] - 33;
        }
    }

    return rp.avgQual;
}

bool isCandidateDup(ReadPair & a, ReadPair & b, unsigned minMatches, unsigned /*maxOffset*/, DupPrefix const)
{
    return prefix(a.read1, minMatches) == prefix(b.read1, minMatches);
}

bool isCandidateDup(ReadPair & a, ReadPair & b, unsigned /*minMatches*/, unsigned maxOffset, DupReadName const)
{
    unsigned sizeA = a.qnameSplit.size();
    unsigned sizeB = b.qnameSplit.size();

    // Check if minimum requirements on read name composition are met.
    if (sizeA < 3 || sizeB < 3 || sizeA != sizeB)
        return false;

    // Check if reads are from the same tile.
    if (a.qnameSplit[sizeA-3] != b.qnameSplit[sizeB-3])
        return false;

    // Check if reads are nearby (within maxOffset) on the same tile.
    int ax = atoi(seqan::toCString(a.qnameSplit[sizeA-2]));
    int bx = atoi(seqan::toCString(b.qnameSplit[sizeB-2]));
    if (abs(ax - bx) < maxOffset)
        return true;
    else
    {
        int ay = atoi(seqan::toCString(a.qnameSplit[sizeA-1]));
        int by = atoi(seqan::toCString(b.qnameSplit[sizeB-1]));
        return abs(ay - by) < maxOffset;
    }
}

bool isCandidateDup(ReadPair & a, ReadPair & b, unsigned minMatches, unsigned maxOffset, DupBoth const)
{
    return isCandidateDup(a, b, minMatches, maxOffset, DupPrefix()) ||
        isCandidateDup(a, b, minMatches, maxOffset, DupReadName());
}