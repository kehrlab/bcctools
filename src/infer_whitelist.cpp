#include <fstream>
#include <seqan/stream.h>
#include "infer_whitelist.h"
#include "utils.h"

using namespace seqan;

void make_histograms(std::vector<unsigned> & allHist, std::vector<unsigned> & wlHist, std::vector<uint16_t> & count_per_barcode, CharString & whitelistFile, unsigned bcLength)
{
    resize(allHist, 1000, 0u);
    resize(wlHist, 1000, 0u);
    
    // Read the whitelist file.
    std::vector<bool> whitelisted;
    resize(whitelisted, (uint64_t)1 << (2*bcLength), 0);
    
    if (whitelistFile != "")
    {
        printStatus("Reading whitelist file");
        std::ifstream infile(toCString(whitelistFile));
        std::string barcode;
        DnaString bc;
        while (infile >> barcode)
        {
            bc = barcode;
            whitelisted[hash(bc)] = 1;
        }
        printDone();
    }
    
    // Count the barcodes in histograms.
    printStatus("Making barcode histogram");
    for (uint64_t i = 0; i < count_per_barcode.size(); ++i)
    {
        unsigned cnt = 1000;
        if (count_per_barcode[i] < cnt)
            cnt = count_per_barcode[i];
        ++allHist[cnt];
        if (whitelistFile != "" && whitelisted[i])
            ++wlHist[cnt];
    }
    printDone();
}

double entropy(DnaString & bc)
{
    // Count dinucleotide occurrences
    String<unsigned> diCounts;
    resize(diCounts, 16, 0);
    for (unsigned i = 0; i < length(bc)-1; ++i)
    {
        diCounts[ordValue(bc[i]) + 4*ordValue(bc[i+1])] += 1;
    }

    // Calculate entropy score for dinucleotide counts
    double score = 0.0;
    for (unsigned i = 0; i < length(diCounts); ++i)
    {
        if (diCounts[i] == 0) continue;
        double p = double(diCounts[i]) / (length(bc) - 1);
        score -= p * log(p) / log(2);
    }

    return score / 4;
}

unsigned infer_cutoff(std::vector<unsigned> & allHist, std::vector<unsigned> & /*wlHist*/)
{
    unsigned cutoff = 1;
    
    unsigned min = maxValue<unsigned>()/2;
    for (unsigned i = 0; i < allHist.size(); ++i)
    {
        if (allHist[i] < min)
        {
            min = allHist[i];
            cutoff = i;
        }
        if (allHist[i] > 2*min)
            break;
    }
    
    return cutoff;
}
