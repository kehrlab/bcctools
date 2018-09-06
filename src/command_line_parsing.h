#ifndef COMMAND_LINE_PARSING_H_
#define COMMAND_LINE_PARSING_H_

#include <seqan/arg_parse.h>

enum Command
{
    BC_INFER_WHITELIST = 0,
    BC_INDEX = 1,
    BC_CORRECT = 2,
    BC_STATS = 3
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

    Options() :
        bcLength(16), spacerLength(7), whitelistCutoff(0), minEntropy(0.5), numAlts(4)
    {}
};

bool fileExists(seqan::CharString const & filename);
seqan::ArgumentParser::ParseResult parseCommandLine(Options & options, int argc, char const ** argv);

#endif  // COMMAND_LINE_PARSING_H_
