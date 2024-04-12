#include <string>
#include <seqan/arg_parse.h>

#include "command_line_parsing.h"
#include "utils.h"

using namespace seqan;

// Returns true if file exists, otherwise false.
bool fileExists(CharString const & filename)
{
    struct stat buffer;
    return (stat(toCString(filename), &buffer) == 0);
}

// Returns true if directory exists, otherwise false.
bool dirExists(CharString const & filename)
{
    std::string str = toCString(filename);
    size_t pos = str.find_last_of("/\\");
    if (pos != std::string::npos)
    {
        struct stat buffer;
        return (stat(toCString(str.substr(0, pos)), &buffer) == 0 && buffer.st_mode && S_IFDIR);
    }
    return true;
}

void addOptionsWhitelist(ArgumentParser & parser, Options & options)
{
    addOption(parser, ArgParseOption("c", "cutoff", "Minimum number of occurences for including barcode in whitelist. Default: Inferred.", ArgParseArgument::INTEGER));
    addOption(parser, ArgParseOption("e", "entropy", "Minimum entropy of inferred barcode.", ArgParseArgument::DOUBLE));
    setDefaultValue(parser, "entropy", options.minEntropy);
    setMinValue(parser, "entropy", "0.0");
    setMaxValue(parser, "entropy", "1.0");

    addOption(parser, ArgParseOption("w", "whitelist", "Whitelist file for making a detailed barcode counts histogram.", ArgParseArgument::INPUT_FILE));
    addOption(parser, ArgParseOption("o", "out", "Name of whitelist output file.", ArgParseArgument::OUTPUT_FILE));
    setDefaultValue(parser, "out", "barcode_whitelist.txt");

    addOption(parser, ArgParseOption("b", "bc-len", "Length of barcode sequence.", ArgParseArgument::INTEGER));
    setDefaultValue(parser, "bc-len", options.bcLength);
    setMinValue(parser, "bc-len", "1");
    setAdvanced(parser, "bc-len");
}

void addAdvancedOptionsWhitelist(ArgumentParser & /*parser*/, Options & /*options*/)
{
    // addOption(parser, ArgParseOption("l", "long", "Description", ArgParseArgument::INTEGER));
    // setAdvanced(parser, "long");
}

void setupParserWhitelist(ArgumentParser & parser, Options & options)
{
    setVersion(parser, VERSION);
    setDate(parser, DATE);
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fIFASTQ1\\fP");

    addDescription(parser, "Infers a whitelist of barcodes from the occurrences of barcodes in the given input FASTQ "
        "file. Barcodes are included in the whitelist if they label more than a fixed number of read pairs. This "
        "number is set to the first local minimum of the barcode occurrence histogram by default or can be specified"
        "with the --cutoff option. Barcodes with an entropy lower than specified with the --entropy option are "
        "excluded from the whitelist but counted in the histogram.");

    // Define the required arguments.
    ArgParseArgument arg1(ArgParseArgument::INPUT_FILE, "FASTQ1", false);
    setHelpText(arg1, "File in FASTQ format containing the first reads in pairs.");
    setValidValues(arg1, "fq fastq FQ FASTQ fq.gz fastq.gz FQ.gz FASTQ.gz");
    addArgument(parser, arg1);

    // Add options and advanced options. The latter are only visible in the full help.
    addOptionsWhitelist(parser, options);
    addAdvancedOptionsWhitelist(parser, options);
}

void addOptionsIndex(ArgumentParser & parser, Options & options)
{
    addOption(parser, ArgParseOption("a", "alts", "Maximum number of alternative corrections.", ArgParseArgument::INTEGER));
    setDefaultValue(parser, "alts", options.numAlts);
    setMinValue(parser, "alts", "1");
    setMaxValue(parser, "alts", "48");
}

void addAdvancedOptionsIndex(ArgumentParser & /*parser*/, Options & /*options*/)
{
    // addOption(parser, ArgParseOption("l", "long", "Description", ArgParseArgument::INTEGER));
    // setAdvanced(parser, "long");
}

void setupParserIndex(ArgumentParser & parser, Options & options)
{
    setVersion(parser, VERSION);
    setDate(parser, DATE);
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fIWHITELIST\\fP");

    addDescription(parser, "Constructs an index of a barcode whitelist file. The --alts option determines the maximum "
        "number of alternative barcode corrections (= whitelisted barcodes) that will be found when using the "
        "ouput index in the 'correct' command. Any barcode that has more alternative corrections with the same "
        "number of substitutions will not be corrected. The value specified with the --alts option will be rounded to "
        "the next larger power of 2. With larger values, index construction takes longer and the index takes up more "
        "space.");

    // Define the required arguments.
    ArgParseArgument arg1(ArgParseArgument::INPUT_FILE, "WHITELIST", false);
    setHelpText(arg1, "File containing barcode whitelist, one barcode per line.");
    addArgument(parser, arg1);

    // Add options and advanced options. The latter are only visible in the full help.
    addOptionsIndex(parser, options);
    addAdvancedOptionsIndex(parser, options);
}

void addOptionsCorrect(ArgumentParser & parser, Options & options)
{
    addOption(parser, ArgParseOption("a", "alts", "Maximum number of alternative corrections.", ArgParseArgument::INTEGER));
    setDefaultValue(parser, "alts", options.numAlts);
    setMinValue(parser, "alts", "1");
    setMaxValue(parser, "alts", "48");

    addOption(parser, ArgParseOption("s", "spacer", "Length of spacer between barcode and read sequence.", ArgParseArgument::INTEGER));
    setDefaultValue(parser, "spacer", options.spacerLength);
}

void addAdvancedOptionsCorrect(ArgumentParser & /*parser*/, Options & /*options*/)
{
    // addOption(parser, ArgParseOption("l", "long", "Description", ArgParseArgument::INTEGER));
    // setAdvanced(parser, "long");
}

void setupParserCorrect(ArgumentParser & parser, Options & options)
{
    setVersion(parser, VERSION);
    setDate(parser, DATE);
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fIWHITELIST\\fP \\fIFASTQ1\\fP \\fIFASTQ2\\fP");

    addDescription(parser, "Extracts barcodes from reads and corrects them using an index of a barcode whitelist. "
        "An index of the barcode whitelist is created for correction on-the-fly unless it has been precomputed with "
        "the index command. Correction is possible for barcodes with a Hamming distance of one to a whitelisted "
        "barcode and a maximum number of alternative corrections (= whitelisted barcodes). See the help message for "
        " the 'index' command for further details on the --alts option.");

    // Define the required arguments.
    ArgParseArgument arg1(ArgParseArgument::INPUT_FILE, "WHITELIST", false);
    setHelpText(arg1, "File containing barcode whitelist.");
    addArgument(parser, arg1);
    ArgParseArgument arg2(ArgParseArgument::INPUT_FILE, "FASTQ1", false);
    setHelpText(arg2, "File in FASTQ format containing the first reads in pairs.");
    setValidValues(arg2, "fq fastq FQ FASTQ fq.gz fastq.gz FQ.gz FASTQ.gz");
    addArgument(parser, arg2);
    ArgParseArgument arg3(ArgParseArgument::INPUT_FILE, "FASTQ2", false);
    setHelpText(arg3, "File in FASTQ format containing the second reads in pairs.");
    setValidValues(arg3, "fq fastq FQ FASTQ fq.gz fastq.gz FQ.gz FASTQ.gz");
    addArgument(parser, arg3);

    // Add options and advanced options. The latter are only visible in the full help.
    addOptionsCorrect(parser, options);
    addAdvancedOptionsCorrect(parser, options);
}

void addOptionsStats(ArgumentParser & parser, Options & options)
{
    addOption(parser, ArgParseOption("o", "out", "Filename of output file.", ArgParseArgument::OUTPUT_FILE));
    options.outFile = "stats.txt";
    setDefaultValue(parser, "out", options.outFile);
}

void addAdvancedOptionsStats(ArgumentParser & /*parser*/, Options & /*options*/)
{
    // addOption(parser, ArgParseOption("l", "long", "Description", ArgParseArgument::INTEGER));
    // setAdvanced(parser, "long");
}

void setupParserStats(ArgumentParser & parser, Options & options)
{
    setVersion(parser, VERSION);
    setDate(parser, DATE);
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fIFASTQ1\\fP");

    addDescription(parser, "Computes simple barcode statistics from a bcctools TSV file or a barcoded FASTQ, SAM, or "
        "BAM file.");

    // Define the required arguments.
    ArgParseArgument arg1(ArgParseArgument::INPUT_FILE, "INPUT_FILE", false);
    setHelpText(arg1, "File in TSV/(gzipped) FASTQ/SAM/BAM format containing barcode-corrected read pairs. In case of "
        "FASTQ specify only file of first reads in pair.");
    setValidValues(arg1, "fq fastq fq.gz fastq.gz sam bam tsv");
    addArgument(parser, arg1);

    // Add options and advanced options. The latter are only visible in the full help.
    addOptionsStats(parser, options);
    addAdvancedOptionsStats(parser, options);
}

void addOptionsDedup(ArgumentParser & parser, Options & options)
{
    addOption(parser, ArgParseOption("m", "minMatches", "Minimum number of matching bases at the beginning of two reads in order to be considered as potential duplicates for sequence based duplicate detection.", ArgParseArgument::INTEGER));
    setDefaultValue(parser, "minMatches", options.minMatches);
    setMinValue(parser, "minMatches", "1");

    addOption(parser, ArgParseOption("o", "maxOffset", "Maximum offset of x and y coordinates on same tile for two read names in order to be considered as potential duplicates for read name based duplicate detection.", ArgParseArgument::INTEGER));
    setDefaultValue(parser, "maxOffset", options.maxOffset);
    setMinValue(parser, "maxOffset", "1");

    addOption(parser, ArgParseOption("d", "maxDiffRate", "Maxixmum number of differences per read length for alignment.", ArgParseArgument::DOUBLE));
    setDefaultValue(parser, "maxDiffRate", options.maxDiffRate);
    //setMinValue(parser, "maxDiffRate", "0.001");

    addOption(parser, ArgParseOption("q", "minQual", "Minimum quality value for a position in the quality string to be considered for total read quality.", ArgParseArgument::INTEGER));
    setDefaultValue(parser, "minQual", options.minQual);
    setMinValue(parser, "minQual", "0");
    setMaxValue(parser, "minQual", "42");

    addOption(parser, ArgParseOption("n", "nameDups", "Mark only duplicates with similar coordinates on tile."));

    addOption(parser, ArgParseOption("s", "seqDups", "Mark only duplicates with matching read1 prefix."));
}

void addAdvancedOptionsDedup(ArgumentParser & /*parser*/, Options & /*options*/)
{
    // addOption(parser, ArgParseOption("l", "long", "Description", ArgParseArgument::INTEGER));
    // setAdvanced(parser, "long");
}

void setupParserDedup(ArgumentParser & parser, Options & options)
{
    setVersion(parser, VERSION);
    setDate(parser, DATE);
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fITSV\\fP");

    addDescription(parser, "Marks duplicate reads similar to 'Picard MarkDuplicates'. Takes as input a TSV file as "
        "written by the correct command that needs to be sorted by barcode.");

    // Define the required arguments.
    ArgParseArgument arg1(ArgParseArgument::INPUT_FILE, "INPUT_FILE", false);
    setValidValues(arg1, "tsv");
    addArgument(parser, arg1);

    // Add options and advanced options. The latter are only visible in the full help.
    addOptionsDedup(parser, options);
    addAdvancedOptionsDedup(parser, options);
}

void getArgumentValuesWhitelist(Options & options, ArgumentParser & parser)
{
    getArgumentValue(options.fastqFile1, parser, 0);
}

void getArgumentValuesIndex(Options & options, ArgumentParser & parser)
{
    getArgumentValue(options.whitelistFile, parser, 0);
}

void getArgumentValuesCorrect(Options & options, ArgumentParser & parser)
{
    getArgumentValue(options.whitelistFile, parser, 0);
    getArgumentValue(options.fastqFile1, parser, 1);
    getArgumentValue(options.fastqFile2, parser, 2);
}

void getArgumentValuesStats(Options & options, ArgumentParser & parser)
{
    getArgumentValue(options.inputFile, parser, 0);
}

void getArgumentValuesDedup(Options & options, ArgumentParser & parser)
{
    getArgumentValue(options.inputFile, parser, 0);
}

void getOptionValuesWhitelist(Options & options, ArgumentParser & parser)
{
    getOptionValue(options.whitelistCutoff, parser, "cutoff");
    getOptionValue(options.whitelistFile, parser, "whitelist");
    getOptionValue(options.outFile, parser, "out");
    getOptionValue(options.bcLength, parser, "bc-len");
}

void getOptionValuesIndex(Options & options, ArgumentParser & parser)
{
    getOptionValue(options.numAlts, parser, "alts");
}

void getOptionValuesCorrect(Options & options, ArgumentParser & parser)
{
    getOptionValue(options.numAlts, parser, "alts");
    getOptionValue(options.spacerLength, parser, "spacer");
}

void getOptionValuesStats(Options & options, ArgumentParser & parser)
{
    getOptionValue(options.outFile, parser, "out");
}

void getOptionValuesDedup(Options & options, ArgumentParser & parser)
{
    getOptionValue(options.minMatches, parser, "minMatches");
    getOptionValue(options.maxOffset, parser, "maxOffset");
    getOptionValue(options.maxDiffRate, parser, "maxDiffRate");
    getOptionValue(options.minQual, parser, "minQual");

    if(isSet(parser, "seqDups"))
        options.nameDups = false;

    if(isSet(parser, "nameDups"))
        options.seqDups = false;
}

ArgumentParser::ParseResult checkOptionValuesWhitelist(Options & options)
{
    SEQAN_TRY
    {
        std::stringstream what;

        if (!fileExists(options.fastqFile1))
        {
            what << "The first input FASTQ file '" << options.fastqFile1 << "' does not exist.";
            SEQAN_THROW(ParseError(what.str()));
        }

	if (!dirExists(options.outFile))
        {
            what << "The path to the output file '" << options.outFile << "' does not exist.";
            SEQAN_THROW(ParseError(what.str()));
        }

        if (options.whitelistFile != "" && !fileExists(options.whitelistFile))
        {
            what << "The given barcode whitelist file '" << options.whitelistFile << "' does not exist.";
            SEQAN_THROW(ParseError(what.str()));
        }
    }
    SEQAN_CATCH(ParseError & ex)
    {
        std::cerr << "ERROR: " << ex.what() << std::endl;
        return ArgumentParser::PARSE_ERROR;
    }
    return ArgumentParser::PARSE_OK;
}

ArgumentParser::ParseResult checkOptionValuesIndex(Options & options)
{
    SEQAN_TRY
    {
        std::stringstream what;
        unsigned numAltsBase = std::ceil(std::log(options.numAlts)/std::log(2));
        if (options.numAlts != 1u << numAltsBase)
        {
            options.numAlts = 1u << numAltsBase;
            std::ostringstream msg;
            msg << "Maximum number of alternative corrections raised to " << options.numAlts << ".";
            printWarning(msg);
        }

        if (!fileExists(options.whitelistFile))
        {
            what << "The given barcode whitelist file '" << options.whitelistFile << "' does not exist.";
            SEQAN_THROW(ParseError(what.str()));
        }
    }
    SEQAN_CATCH(ParseError & ex)
    {
        std::cerr << "ERROR: " << ex.what() << std::endl;
        return ArgumentParser::PARSE_ERROR;
    }
    return ArgumentParser::PARSE_OK;
}

ArgumentParser::ParseResult checkOptionValuesCorrect(Options & options)
{
    SEQAN_TRY
    {
        std::stringstream what;
        unsigned numAltsBase = std::ceil(std::log(options.numAlts)/std::log(2));
        if (options.numAlts != 1u << numAltsBase)
        {
            options.numAlts = 1u << numAltsBase;
            std::ostringstream msg;
            msg << "Maximum number of alternative corrections raised to " << options.numAlts << ".";
            printWarning(msg);
        }

        if (options.spacerLength < 0)
        {
            what << "The given spacer length " << options.spacerLength << " is not in the range of allowed values [spacer length >= 0].";
            SEQAN_THROW(ParseError(what.str()));
        }

        if (!fileExists(options.whitelistFile))
        {
            what << "The given barcode whitelist file '" << options.whitelistFile << "' does not exist.";
            SEQAN_THROW(ParseError(what.str()));
        }

        if (!fileExists(options.fastqFile1))
        {
            what << "The first input FASTQ file '" << options.fastqFile1 << "' does not exist.";
            SEQAN_THROW(ParseError(what.str()));
        }

        if (!fileExists(options.fastqFile2))
        {
            what << "The second input FASTQ file '" << options.fastqFile2 << "' does not exist.";
            SEQAN_THROW(ParseError(what.str()));
        }
    }
    SEQAN_CATCH(ParseError & ex)
    {
        std::cerr << "ERROR: " << ex.what() << std::endl;
        return ArgumentParser::PARSE_ERROR;
    }
    return ArgumentParser::PARSE_OK;
}

ArgumentParser::ParseResult checkOptionValuesStats(Options & options)
{
    SEQAN_TRY
    {
        std::stringstream what;
        if (!fileExists(options.inputFile))
        {
            what << "The input file '" << options.inputFile << "' does not exist.";
            SEQAN_THROW(ParseError(what.str()));
        }

        if (!dirExists(options.outFile))
        {
            what << "The directory of the given output file '" << options.outFile << "' does not exist.";
            SEQAN_THROW(ParseError(what.str()));
        }
    }
    SEQAN_CATCH(ParseError & ex)
    {
        std::cerr << "ERROR: " << ex.what() << std::endl;
        return ArgumentParser::PARSE_ERROR;
    }
    return ArgumentParser::PARSE_OK;
}

ArgumentParser::ParseResult checkOptionValuesDedup(Options & options)
{
    SEQAN_TRY
    {
        std::stringstream what;
        if (!fileExists(options.inputFile))
        {
            what << "The input file '" << options.inputFile << "' does not exist.";
            SEQAN_THROW(ParseError(what.str()));
        }

        if (options.seqDups == false && options.nameDups == false)
        {
            what << "Options --seqDups (-s) and --nameDups (-n) are exclusive. Please specify only one or none.";
            SEQAN_THROW(ParseError(what.str()));
        }
    }
    SEQAN_CATCH(ParseError & ex)
    {
        std::cerr << "ERROR: " << ex.what() << std::endl;
        return ArgumentParser::PARSE_ERROR;
    }
    return ArgumentParser::PARSE_OK;
}

void printHeader(ArgumentParser & parser, std::ostringstream & cmd)
{
    std::ostream_iterator<char> out(std::cerr);
    CharString name = getName(parser._toolDoc);
    CharString shortDescription = getShortDescription(parser._toolDoc);

    // Print the header.
    std::cerr << name;
    if (!empty(shortDescription))
        std::cerr << " - " << shortDescription;
    std::cerr << "\n";
    unsigned len = length(name) + (empty(shortDescription) ? 0 : 3) + length(shortDescription);
    std::fill_n(out, len, '=');
    std::cerr << "\n\n";

    // Print the command.
    std::cerr << "Command: \n" << cmd.str() << std::endl;
}

ArgumentParser::ParseResult parseCommandLine(Options & options, int argc, char const ** argv)
{
    // Define function types for function arrays.
    typedef void (*SetupParserFunctionType)(ArgumentParser &, Options &);
    typedef void (*GetValuesFunctionType)(Options &, ArgumentParser &);
    typedef ArgumentParser::ParseResult (*CheckValuesFunctionType)(Options &);

    // Initialize function arrays.
    SetupParserFunctionType setupParser[5] = {&setupParserWhitelist, &setupParserIndex, &setupParserCorrect, &setupParserStats, &setupParserDedup};
    GetValuesFunctionType getArgumentValues[5] = {&getArgumentValuesWhitelist, &getArgumentValuesIndex, &getArgumentValuesCorrect, &getArgumentValuesStats, &getArgumentValuesDedup};
    GetValuesFunctionType getOptionValues[5] = {& getOptionValuesWhitelist, &getOptionValuesIndex, &getOptionValuesCorrect, &getOptionValuesStats, &getOptionValuesDedup};
    CheckValuesFunctionType checkOptionValues[5] = {&checkOptionValuesWhitelist, &checkOptionValuesIndex, &checkOptionValuesCorrect, &checkOptionValuesStats, &checkOptionValuesDedup};

    // Retrieve the command line.
    std::ostringstream command_line;
    for (int i = 0; i < argc; ++i)
        command_line << argv[i] << " ";
    command_line << std::endl;

    // Store the command.
    const char * command = argv[1];
    if (strcmp(command, "whitelist") == 0)
        options.cmd = Command::BC_INFER_WHITELIST;
    else if (strcmp(command, "index") == 0)
        options.cmd = Command::BC_INDEX;
    else if (strcmp(command, "correct") == 0)
        options.cmd = Command::BC_CORRECT;
    else if (strcmp(command, "stats") == 0)
        options.cmd = Command::BC_STATS;
    else if (strcmp(command, "dedup") == 0)
        options.cmd = Command::BC_DEDUP;
    else
    {
        std::cerr << "ERROR: Unknown command '" << command << "'." << std::endl;
        return ArgumentParser::PARSE_ERROR;
    }

    // Concatenate program name and command.
    const char * p = argv[0];
    std::ostringstream name;
    name << p << " " << command;
    CharString prog_name = name.str();

    ++argv;
    --argc;

    // Setup the parser.
    ArgumentParser parser(toCString(prog_name));
    setupParser[options.cmd](parser, options);

    // Parse the command line and write error messages to error stream.
    std::ostringstream errorStream;
    ArgumentParser::ParseResult res = parse(parser, argc, argv, std::cout, errorStream);

    // If error occurred in parsing, now print the error message and return.
    if (res != ArgumentParser::PARSE_OK)
    {
        std::cerr << errorStream.str();
        return res;
    }

    // Print header and command line to std::out.
    printHeader(parser, command_line);

    // Get the required arguments' and parameters' values.
    getArgumentValues[options.cmd](options, parser);
    getOptionValues[options.cmd](options, parser);

    // Check that parameter values are valid, e.g. input files exist.
    return checkOptionValues[options.cmd](options);
}
