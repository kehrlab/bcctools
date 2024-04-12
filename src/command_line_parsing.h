#ifndef COMMAND_LINE_PARSING_H_
#define COMMAND_LINE_PARSING_H_

#include <seqan/arg_parse.h>

enum Command
{
    BC_INFER_WHITELIST = 0,
    BC_INDEX = 1,
    BC_CORRECT = 2,
    BC_STATS = 3,
    BC_DEDUP = 4
};

struct Options
{
    Command cmd;

    seqan::CharString whitelistFile;
    seqan::CharString fastqFile1;
    seqan::CharString fastqFile2;
    seqan::CharString inputFile;
    seqan::CharString outFile;

    int bcLength;
    int spacerLength;
    unsigned whitelistCutoff;
    double minEntropy;
    unsigned numAlts;

    unsigned minMatches;
    unsigned maxOffset;
    double maxDiffRate;
    unsigned minQual;
    bool nameDups;
    bool seqDups;

    Options() :
        bcLength(16), spacerLength(7), whitelistCutoff(0), minEntropy(0.5), numAlts(16),
        minMatches(5), maxOffset(5000), maxDiffRate(0.1), minQual(15), nameDups(true), seqDups(true)
    {}
};

bool fileExists(seqan::CharString const & filename);
seqan::ArgumentParser::ParseResult parseCommandLine(Options & options, int argc, char const ** argv);

#endif  // COMMAND_LINE_PARSING_H_
