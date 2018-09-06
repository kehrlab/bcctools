#include <unistd.h>
#include <zlib.h>
#include <htslib/kseq.h>

#include "utils.h"
#include "command_line_parsing.h"
#include "infer_whitelist.h"
#include "barcode_index.h"
#include "stats.h"

using namespace seqan;

void printHelp(char const * name)
{
    std::cerr << "bctools - Barcode correction, sorting, etc." << std::endl;
    std::cerr << "===========================================" << std::endl;
    std::cerr << std::endl;
    std::cerr << "\033[1mSYNOPSIS\033[0m" << std::endl;
    std::cerr << "    \033[1m" << name << " COMMAND\033[0m [\033[4mOPTIONS\033[0m]" << std::endl;
    std::cerr << std::endl;
    std::cerr << "\033[1mCOMMANDS\033[0m" << std::endl;
    std::cerr << "    \033[1mwhitelist\033[0m Infers a whitelist file from the number of barcode occurences." << std::endl;
    std::cerr << "    \033[1mindex\033[0m     Builds the barcode index from a barcode whitelist and writes it to files." << std::endl;
    std::cerr << "    \033[1mcorrect\033[0m   Cuts off the barcodes in a pair of FASTQ files and corrects them using a barcode whitelist." << std::endl;
    std::cerr << "    \033[1mstats\033[0m     Computes barcode statistics for the given TSV/FASTQ/SAM/BAM file." << std::endl;
    std::cerr << std::endl;
    std::cerr << "\033[1mVERSION\033[0m" << std::endl;
    std::cerr << "    " << name << " version: " << VERSION << std::endl;
    std::cerr << "    Last update " << DATE << std::endl;
    std::cerr << std::endl;
    std::cerr << "Try `" << name << " COMMAND --help' for more information on each command." << std::endl;
    std::cerr << std::endl;
}

void count_corrected_pair(BarcodeStatus s, unsigned & match, unsigned & one_error, unsigned & unrecognized, unsigned & invalid)
{
    switch (s)
    {
        case BarcodeStatus::MATCH:
            ++match; return;
        case BarcodeStatus::ONE_ERROR:
            ++one_error; return;
        case BarcodeStatus::UNRECOGNIZED:
            ++unrecognized; return;
        case BarcodeStatus::INVALID:
            ++invalid; return;
    }
}

void print_correction_stats(unsigned match, unsigned one_error, unsigned unrecognized, unsigned invalid)
{
    std::cerr << "Stats:" << std::endl;
    std::cerr << "  Whitelisted barcodes:    " << match << std::endl;
    std::cerr << "  Corrected barcodes:  " << one_error << std::endl;
    std::cerr << "  Unrecognized barcodes:   " << unrecognized << std::endl;
    if (invalid > 0)
        std::cerr << "  Barcodes invalid:       " << invalid << std::endl;
}

void write_tsv(kseq_t * seq1, kseq_t * seq2, std::vector<seqan::DnaString> & barcodeCorrected, unsigned bcLength, int spacerLength)
{
    // Field 1: Read name (ID).
    std::cout << seq1->name.s;

    // Field 2: Corrected barcode or '*' if none.
    if (barcodeCorrected.size() != 0)
    {
        std::cout << "\t" << barcodeCorrected[0];
        for (unsigned i = 1; i < barcodeCorrected.size(); ++i)
            std::cout << "," << barcodeCorrected[i];
    }
    else
        std::cout << "\t" << "*";

    // Field 3 and 4: Raw barcode and spacer sequence.
    std::cout << "\t";
    for (unsigned i = 0; i < bcLength; ++i)
        std::cout << seq1->seq.s[i];
    std::cout << "\t";
    for (unsigned i = bcLength; i < bcLength + spacerLength; ++i)
        std::cout << seq1->seq.s[i];

    // Field 5 and 6: Sequence of first and second read.
    std::cout << "\t";
    for (unsigned i = bcLength + spacerLength; i < strlen(seq1->seq.s); ++i)
        std::cout << seq1->seq.s[i];
    std::cout << "\t" << seq2->seq.s;

    // Field 7 and 8: Barcode and spacer qual.
    std::cout << "\t";
    for (unsigned i = 0; i < bcLength; ++i)
        std::cout << seq1->qual.s[i];
    std::cout << "\t";
    for (unsigned i = bcLength; i < bcLength + spacerLength; ++i)
        std::cout << seq1->qual.s[i];

    // Field 9 and 10: Quality string of first and second read.
    std::cout << "\t";
    for (unsigned i = bcLength + spacerLength; i < strlen(seq1->qual.s); ++i)
        std::cout << seq1->qual.s[i];
    std::cout << "\t" << seq2->qual.s;

    std::cout << std::endl;
}

int infer_whitelist(Options & options)
{
    // Open the input FASTQ file.
    printStatus("Opening FASTQ file");
    gzFile fp1 = gzopen(toCString(options.fastqFile1), "r");
    kseq_t * seq1 = kseq_init(fp1);
    printDone();

    // Initialize counts table.
    std::vector<uint16_t> count_per_barcode;
    resize(count_per_barcode, (uint64_t)1 << (2*options.bcLength), 0);

    // Iterate the FASTQ records and count barcodes.
    printStatus("Counting barcodes");
    while (kseq_read(seq1) >= 0)
    {
        Dna5String rx = prefix(seq1->seq.s, options.bcLength);
        typename Iterator<Dna5String, Rooted>::Type it = begin(rx);
        typename Iterator<Dna5String, Rooted>::Type itEnd = end(rx);
        bool hasN = false;
        for (; it != itEnd; ++it)
        {
            if (*it == 'N')
            {
                hasN = true;
                break;
            }
        }
        if (hasN == false)
        {
            DnaString rrx = rx;
            uint64_t h = hash(rrx);
            if (count_per_barcode[h] != maxValue<uint16_t>())
                ++count_per_barcode[h];
        }
    }
    printDone();

    // Cleanup and close all files.
    kseq_destroy(seq1);
    gzclose(fp1);

    // Make histogram of all barcode counts and, optionally, of whitelisted barcode counts.
    std::vector<unsigned> allHist;
    std::vector<unsigned> wlHist;
    make_histograms(allHist, wlHist, count_per_barcode, options.whitelistFile, options.bcLength);

    // Write the histograms to file.
    printStatus("Writing histograms of barcodes");
    CharString filename = options.outFile;
    filename += ".hist";
    std::ofstream histFile(toCString(filename));
    histFile << "All" << "\t" << "Whitelisted" << std::endl;
    for (uint64_t i = 0; i < allHist.size(); ++i)
    {
        histFile << allHist[i] << "\t" << wlHist[i] << std::endl;
    }
    histFile.close();
    printDone();

    // Infer cutoff for whitelisting a barcode.
    if (options.whitelistCutoff == 0)
        options.whitelistCutoff = infer_cutoff(allHist, wlHist);

    std::ostringstream msg;
    msg << "Minimum number of barcode occurences set to '" << options.whitelistCutoff << "'.";
    printInfo(msg);

    // Write the whitelist to file.
    printStatus("Writing whitelist of barcodes");
    std::ofstream outFile(toCString(options.outFile));
    for (uint64_t i = 0; i < count_per_barcode.size(); ++i)
    {
        if (count_per_barcode[i] >= options.whitelistCutoff)
        {
            DnaString bc = unhash(i, options.bcLength);
            double e = entropy(bc);
            if (e >= options.minEntropy)
                outFile << bc << std::endl;
        }
    }
    outFile.close();
    printDone();

    return 0;
}

int index(Options & options)
{
    BarcodeIndex sbi(options.whitelistFile);

    // Build the barcode index.
    buildIndex(sbi, options.whitelistFile, options.numAlts);

    // Write the index to files.
    writeBarcodeTable(options.whitelistFile, sbi);
    writeMatchTable(options.whitelistFile, sbi);
    writeSubstitutionTable(options.whitelistFile, sbi);

    return 0;
}

int correct(Options & options)
{
    BarcodeIndex sbi(options.whitelistFile);

    CharString bcFilename = options.whitelistFile;
    bcFilename += ".bc";
    if (!fileExists(bcFilename))
    {
        // Build the barcode index.
        buildIndex(sbi, options.whitelistFile, options.numAlts);
    }
    else
    {
        // Load the barcode index.
        if (load(sbi, options.whitelistFile) != 0)
        {
            std::ostringstream what;
            what << "Some barcode index files for " << options.whitelistFile << "' do not exist. Use the 'index' command to create them.";
            SEQAN_THROW(ParseError(what.str()));
        }
        std::ostringstream msg;
        msg << "Maximum number of alternative corrections stored in index is " << sbi.numAlts << ".";
        printInfo(msg);
    }

    // Open the input and output files.
    printStatus("Opening FASTQ files");
    gzFile fp1 = gzopen(toCString(options.fastqFile1), "r");
    gzFile fp2 = gzopen(toCString(options.fastqFile2), "r");
    kseq_t * seq1 = kseq_init(fp1);
    kseq_t * seq2 = kseq_init(fp2);
    printDone();

    // Initialize counts on barcode correction.
    unsigned match, one_error, unrecognized, invalid;
    match = one_error = unrecognized = invalid = 0;

    // Iterate the FASTQ records and retrieve the corrected barcodes.
    printInfo("Retrieving whitelist barcodes.");
    while (kseq_read(seq1) >= 0 && kseq_read(seq2) >=0)
    {
        Dna5String rx = prefix(seq1->seq.s, sbi.bcLength);
        CharString qx = prefix(seq1->qual.s, sbi.bcLength);
        std::vector<seqan::DnaString> barcodeCorrected;
        BarcodeStatus s = retrieve(barcodeCorrected, sbi, rx, qx);
        count_corrected_pair(s, match, one_error, unrecognized, invalid);
        write_tsv(seq1, seq2, barcodeCorrected, sbi.bcLength, options.spacerLength);
    }

    // Cleanup and close all files.
    kseq_destroy(seq1);
    kseq_destroy(seq2);
    gzclose(fp1);
    gzclose(fp2);

    // Print counts on barcode correction.
    std::cerr << std::endl;
    print_correction_stats(match, one_error, unrecognized, invalid);
    std::cerr << std::endl;

    return 0;
}

int stats(Options & options)
{
    typedef ModifiedString<CharString, ModView<FunctorLowcase<char> > > TLowcase;
    TLowcase lowcaseFilename(options.inputFile);

    // Initialize the BarcodeStats.
    BarcodeStats stats;

    // Stream over and count each read pair in input file.
    if (suffix(lowcaseFilename, length(lowcaseFilename) - 3) == "tsv")
        stats_tsv(stats, options.inputFile);
    else if (suffix(lowcaseFilename, length(lowcaseFilename) - 5) == "fastq" ||
             suffix(lowcaseFilename, length(lowcaseFilename) - 8) == "fastq.gz")
        stats_fastq(stats, options.inputFile);
    else if (suffix(lowcaseFilename, length(lowcaseFilename) - 3) == "sam" ||
             suffix(lowcaseFilename, length(lowcaseFilename) - 3) == "bam")
        stats_bam(stats, options.inputFile);

    // Write output.
    printStatus("Writing stats to output file");
    write_stats(options.outFile, stats);
    printDone();

    return 0;
}

// =============================================================================
// Function main()
// =============================================================================

int main(int argc, const char ** argv)
{
    // Turn off synchronization of cin with stdio.
    std::ios_base::sync_with_stdio(false);

    // Check if simply the help message is requested.
    if (argc < 2 || strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0)
    {
        const char * prog_name = argv[0];
        printHelp(prog_name);
        return 1;
    }

    try
    {
        // Parse the command line.
        Options options;
        ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
        if (res == ArgumentParser::PARSE_HELP || res == ArgumentParser::PARSE_VERSION ||
            res ==  ArgumentParser::PARSE_WRITE_CTD || res == ArgumentParser::PARSE_EXPORT_HELP)
            return 0;
        else if (res != ArgumentParser::PARSE_OK)
            return 1;

        // Execute the specified command.
        bool ret = 1;
        if (options.cmd == Command::BC_INFER_WHITELIST)
            ret = infer_whitelist(options);
        else if (options.cmd == Command::BC_INDEX)
            ret = index(options);
        else if (options.cmd == Command::BC_CORRECT)
            ret = correct(options);
        else if (options.cmd == Command::BC_STATS)
            ret = stats(options);

        if (ret == 0)
            printInfo("Finished successfully.");
        return ret;
    }
    catch (seqan::Exception const & ex)
    {
        std::cerr << "ERROR: " << ex.what() << std::endl;
        return 1;
    }
}
