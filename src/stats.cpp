#include <fstream>
#include <seqan/stream.h>
#include <seqan/bam_io.h>
#include "stats.h"
#include "utils.h"

using namespace seqan;

void count_read_pair(BarcodeStats & stats, Dna5String barcode, std::vector<DnaString> barcodeCorrected, const char * qual)
{
    // Count read pair with unrecognized barcode.
    if (barcodeCorrected.size() == 0)
    {
        ++stats.unrecognized;
        return;
    }
    // If not already done, initialize barcode counts histogram.
    stats.count_hist.resize(1000);

    // Count corrected barcode in histogram of barcode occurences.
    if (barcodeCorrected[0] == stats.prev_corrected)
    {
        ++stats.prev_count;
    }
    else
    {
        if (stats.prev_count != 0)
        {
            if (stats.prev_count < stats.count_hist.size() - 1)
                ++stats.count_hist[stats.prev_count];
            else
                ++stats.count_hist[stats.count_hist.size() - 1];
        }
        stats.prev_corrected = barcodeCorrected[0];
        stats.prev_count = 1;
    }

    // Find mismatch positions if any.
    std::vector<unsigned> mismatch_pos;
    for (unsigned i = 0; i < length(barcode); ++i)
    {
        if ((Dna5)barcodeCorrected[0][i] != barcode[i])
            mismatch_pos.push_back(i);
    }

    // Count mismatch positions.
    if (mismatch_pos.size() == 0)
    {
        ++stats.error_free;
    }
    else if (mismatch_pos.size() == 1)
    {
        ++stats.one_error;
        char qv = qual[mismatch_pos[0]];
        if (stats.one_error_hist.count(qv) == 0)
            stats.one_error_hist[qv].resize(length(barcode));
        ++stats.one_error_hist[qv][mismatch_pos[0]];
    }
}

void write_stats(CharString & outputFile, BarcodeStats & stats)
{
    std::ofstream out(toCString(outputFile));

    out << "ERROR_FREE_BARCODES" << "\t" << stats.error_free << std::endl;
    out << "ONE_MISMATCH_BARCODES" << "\t" << stats.one_error << std::endl;
    out << "UNRECOGNIZED_BARCODES" << "\t" << stats.unrecognized << std::endl;

    out << "BARCODE_COUNT_HIST";
    for (unsigned i = 0; i < stats.count_hist.size(); ++i)
        out << "\t" << stats.count_hist[i];
    out << std::endl;

    for (std::map<char, std::vector<unsigned> >::iterator it = stats.one_error_hist.begin(); it != stats.one_error_hist.end(); ++it)
    {
        out << "ONE_ERROR_HIST_QUAL_" << it->first;
        for (unsigned i = 0; i < it->second.size(); ++i)
            out << "\t" << it->second[i];
        out << std::endl;
    }

    out.close();
}

void parse_kseq_comment(Dna5String & barcode, std::vector<seqan::DnaString> & barcodeCorrected, kstring_t * comment)
{
    for (unsigned i = 0; i < strlen(comment->s) - 5; ++i)
    {
        if (comment->s[i] == 'B' && comment->s[i+1] == 'X')
        {
            i += 5;
            if (comment->s[i] == '*')
            {
                ++i;
            }
            else
            {
                while (i < strlen(comment->s) && !isspace(comment->s[i]))
                {
                    DnaString bc;
                    while (i < strlen(comment->s) && !isspace(comment->s[i]) && comment->s[i] != ',')
                    {
                        appendValue(bc, comment->s[i]);
                        ++i;
                    }
                    barcodeCorrected.push_back(bc);
                    if (comment->s[i] == ',')
                        ++i;
                }
            }
        }
        else if (comment->s[i] == 'R' && comment->s[i+1] == 'X')
        {
            i += 5;
            while (i < strlen(comment->s) && !isspace(comment->s[i]))
            {
                appendValue(barcode, comment->s[i]);
                ++i;
            }
        }
    }
}

std::vector<DnaString> parseBarcodeList(const char * cBarcode)
{
    std::vector<DnaString> barcodes;

    if (cBarcode[0] == '*')
        return barcodes;

    std::string bc;
    std::istringstream barcodeStream(cBarcode);
    while (std::getline(barcodeStream, bc, ','))
    {
        DnaString barcode = bc;
        barcodes.push_back(barcode);
    }
    return barcodes;
}

void stats_tsv(BarcodeStats & stats, CharString & inputFile)
{
    std::ifstream in(toCString(inputFile));

    printStatus("Streaming over the input TSV file");

    std::string qname, cBarcode, rBarcode, spacer, read1, read2, qBarcode, qSpacer, qual1, qual2;
    while(in >> qname >> cBarcode >> rBarcode >> spacer >> read1 >> read2 >> qBarcode >> qSpacer >> qual1 >> qual2)
    {
        Dna5String barcode = rBarcode;
        std::vector<DnaString> barcodeCorrected;
        if (length(cBarcode) > 0)
            barcodeCorrected = parseBarcodeList(toCString(cBarcode));
        count_read_pair(stats, rBarcode, barcodeCorrected, toCString(qBarcode));
    }
    in.close();

    printDone();
}

void stats_fastq(BarcodeStats & stats, CharString & inputFile)
{
    // Open input FASTQ file (only first read in pair needed).
    gzFile fp1 = gzopen(toCString(inputFile), "r");
    kseq_t * seq = kseq_init(fp1);

    if (kseq_read(seq) < 0)
        SEQAN_THROW(ParseError("Input FASTQ file is empty."));

    printStatus("Streaming over the input FASTQ file");

    Dna5String bc;
    std::vector<seqan::DnaString> bcCorrected;
    parse_kseq_comment(bc, bcCorrected, &seq->comment);

    // Count the first pair.
    count_read_pair(stats, bc, bcCorrected, seq->qual.s);

    // Iterate and count the remaining FASTQ records.
    while (kseq_read(seq) >= 0)
    {
        Dna5String barcode;
        std::vector<seqan::DnaString> barcodeCorrected;
        parse_kseq_comment(barcode, barcodeCorrected, &seq->comment);
        count_read_pair(stats, barcode, barcodeCorrected, seq->qual.s);
    }

    // Close input file.
    kseq_destroy(seq);
    gzclose(fp1);

    printDone();
}

void stats_bam(BarcodeStats & stats, CharString & inputFile)
{
    BamFileIn in(toCString(inputFile));
    BamHeader header;
    readHeader(header, in);

    printStatus("Streaming over the input SAM/BAM file");

    BamAlignmentRecord record;
    while (!atEnd(in))
    {
        readRecord(record, in);
        
        if (!hasFlagFirst(record))
            continue;
        
        BamTagsDict tagsDict(record.tags);
        unsigned tagIdx = 0;

        CharString cBarcode;
        if (findTagKey(tagIdx, tagsDict, "BX"))
            extractTagValue(cBarcode, tagsDict, tagIdx);
        std::vector<DnaString> barcodeCorrected;
        if (length(cBarcode) > 0)
            barcodeCorrected = parseBarcodeList(toCString(cBarcode));

        Dna5String barcode;
        if (findTagKey(tagIdx, tagsDict, "RX"))
            extractTagValue(barcode, tagsDict, tagIdx);

        CharString qBarcode;
        if (findTagKey(tagIdx, tagsDict, "QX"))
            extractTagValue(qBarcode, tagsDict, tagIdx);

        count_read_pair(stats, barcode, barcodeCorrected, toCString(qBarcode));
    }

    printDone();
}
