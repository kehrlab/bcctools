#include <sdsl/bit_vectors.hpp>
#include <seqan/stream.h>
#include <seqan/sequence.h>

#include "barcode_index.h"
#include "utils.h"

using namespace seqan;

inline BarcodeStatus get_status(BarcodeIndex & sbi, uint64_t h)
{
    if (sbi.barcode_table[h])
    {
        if (sbi.match_table[sbi.rank_support_barcode_table(h)])
            return BarcodeStatus::ONE_ERROR;
        else
            return BarcodeStatus::MATCH;
    }
    else
    {
        return BarcodeStatus::UNRECOGNIZED;
    }
}

inline BarcodeStatus get_status_uncondensed(BarcodeIndex & sbi, uint64_t h)
{
    if (sbi.barcode_table[h])
    {
        if (sbi.match_table[h])
            return BarcodeStatus::ONE_ERROR;
        else
            return BarcodeStatus::MATCH;
    }
    else
    {
        if (sbi.match_table[h])
            return BarcodeStatus::INVALID;
        else
            return BarcodeStatus::UNRECOGNIZED;
    }
}

inline void set_match(BarcodeIndex & sbi, uint64_t h)
{
    sbi.barcode_table[h] = 1;
    sbi.match_table[h] = 0;
}

inline void set_one_error(BarcodeIndex & sbi, uint64_t h, sdsl::int_vector<> & helper_table)
{
    BarcodeStatus current = get_status_uncondensed(sbi, h);
    if (current == BarcodeStatus::UNRECOGNIZED && helper_table[h] != 1)
    {
        // Set to ONE_ERROR.
        sbi.barcode_table[h] = 1;
        sbi.match_table[h] = 1;
        helper_table[h] = 0;
    }
    else if (current == BarcodeStatus::ONE_ERROR && helper_table[h] != sbi.numAlts-1)
    {
        ++helper_table[h];
    }
    else if (current == BarcodeStatus::ONE_ERROR && helper_table[h] == sbi.numAlts-1)
    {
        // Set to INVALID.
        sbi.barcode_table[h] = 0;
        sbi.match_table[h] = 0;
        helper_table[h] = 1;
        // INVALID and helper_table[h] = 1 indicates that too many one_error corrections are possible.
    }
}

inline void add_similar_barcodes(BarcodeIndex & sbi, uint64_t h, sdsl::int_vector<> & helper_table)
{
    for (uint64_t i = 0; i < sbi.bcLength; ++i)
    {
        h ^= static_cast<uint64_t>(1) << 2*i;
        set_one_error(sbi, h, helper_table);

        h ^= static_cast<uint64_t>(2) << 2*i;
        set_one_error(sbi, h, helper_table);

        h ^= static_cast<uint64_t>(1) << 2*i;
        set_one_error(sbi, h, helper_table);

        h ^= static_cast<uint64_t>(2) << 2*i;
    }
}

inline bool get_substitution(unsigned & i, BarcodeIndex & sbi, uint64_t h, unsigned offset)
{
    uint64_t index = (sbi.rank_support_match_table(sbi.rank_support_barcode_table(h)) << (sbi.numAltsBase)) + offset;

    if (offset > 0 && sbi.substitution_table[index] == sbi.substitution_table[index - 1])
        return false;

    i = sbi.substitution_table[index];
    return true;
}


inline bool set_substitution(BarcodeIndex & sbi, uint64_t h, uint64_t i, sdsl::int_vector<> & helper_subst_table)
{
    if (get_status(sbi, h) != BarcodeStatus::ONE_ERROR)
        return 1;

    uint64_t pos = sbi.rank_support_match_table(sbi.rank_support_barcode_table(h));
    unsigned offset = helper_subst_table[pos];
    uint64_t index = (pos << (sbi.numAltsBase)) + offset;

    sbi.substitution_table[index] = i;
    ++helper_subst_table[pos];
    if (offset < sbi.numAlts - 1)
    {
        sbi.substitution_table[index+1] = i;
    }
    return 0;
}

inline uint64_t get_corrected_barcode(BarcodeIndex & sbi, uint64_t h, unsigned i)
{
    h ^= static_cast<uint64_t>(1) << 2*i;
    if (get_status(sbi, h) == BarcodeStatus::MATCH)
        return h;

    h ^= static_cast<uint64_t>(2) << 2*i;
    if (get_status(sbi, h) == BarcodeStatus::MATCH)
        return h;

    h ^= static_cast<uint64_t>(1) << 2*i;
    if (get_status(sbi, h) == BarcodeStatus::MATCH)
        return h;

    h ^= static_cast<uint64_t>(2) << 2*i;
    return h;
}

BarcodeIndex::BarcodeIndex(CharString & filename)
{
    // Open barcode whitelist file and determine barcode length.
    std::ifstream infile(toCString(filename));
    std::string barcode;
    if (!(infile >> barcode))
    {
        std::stringstream what;
        what << "Barcode whitelist file '" << toCString(filename) << "' is empty.";
        SEQAN_THROW(ParseError(what.str()));
    }
    bcLength = barcode.size();
}

void condenseMatchTable(BarcodeIndex & sbi)
{
    printStatus("Condensing match table");
    uint64_t m = 0;
    for (uint64_t b = 0; b < sbi.barcode_table.size(); ++b)
    {
        if (sbi.barcode_table[b] == 1)
        {
            sbi.match_table[m] = sbi.match_table[b];
            ++m;
        }
    }
    sbi.match_table.resize(m);
    printDone();
}

void buildBarcodeAndMatchTable(BarcodeIndex & sbi, seqan::CharString & filename)
{
    printStatus("Building barcode table and match table of barcode index");

    // Initialize the barcode table.
    sbi.barcode_table = sdsl::bit_vector((uint64_t)1 << (2*sbi.bcLength), 0u);
    sbi.match_table = sdsl::bit_vector((uint64_t)1 << (2*sbi.bcLength), 0u);
    sdsl::int_vector<> helper_table((uint64_t)1 << (2*sbi.bcLength), 0u, sbi.numAlts);

    // Fill the barcode table.
    std::ifstream infile(toCString(filename));
    std::string barcode;
    DnaString bc;
    uint64_t h;
    while (infile >> barcode)
    {
        bc = barcode;
        h = hash(bc);
        set_match(sbi, h);
        add_similar_barcodes(sbi, h, helper_table);
    }
    printDone();

    condenseMatchTable(sbi);
}

void buildSubstitutionTable(BarcodeIndex & sbi, seqan::CharString & filename)
{
    printStatus("Building substitution table of barcode index");

    // Initialize the rank support tables.
    sbi.rank_support_barcode_table = sdsl::rank_support_v<>(&sbi.barcode_table);
    sbi.rank_support_match_table = sdsl::rank_support_v<>(&sbi.match_table);

    // Initialize the substitution table.
    unsigned bcPos = std::ceil(std::log(sbi.bcLength) / std::log(2));
    uint64_t subst = sbi.rank_support_match_table(sbi.match_table.size());
    sbi.substitution_table = sdsl::int_vector<>(subst * sbi.numAlts, 0u, bcPos);
    sdsl::int_vector<> helper_subst_table(subst, 0u, sbi.numAlts);

    // Fill the substitution table.
    std::ifstream infile(toCString(filename));
    std::string barcode;
    while (infile >> barcode)
    {
        DnaString bc = barcode;
        uint64_t h = hash(bc);
        for (unsigned i = 0; i < sbi.bcLength; ++i)
        {
            h ^= static_cast<uint64_t>(1) << 2*i;
            set_substitution(sbi, h, i, helper_subst_table);
            h ^= static_cast<uint64_t>(2) << 2*i;
            set_substitution(sbi, h, i, helper_subst_table);
            h ^= static_cast<uint64_t>(1) << 2*i;
            set_substitution(sbi, h, i, helper_subst_table);
            h ^= static_cast<uint64_t>(2) << 2*i;
        }
    }
    printDone();
}

void buildIndex(BarcodeIndex & sbi, seqan::CharString & filename, unsigned alternatives)
{
    sbi.numAlts = alternatives;
    sbi.numAltsBase = std::log(sbi.numAlts) / std::log(2);

    buildBarcodeAndMatchTable(sbi, filename);
    buildSubstitutionTable(sbi, filename);
}

int load(BarcodeIndex & sbi, seqan::CharString & filename)
{
    printStatus("Loading barcode index");

    CharString bcFilename = filename;
    bcFilename += ".bc";

    // Open the barcode table file.
    std::ifstream in;
    in.open(toCString(bcFilename), std::ios::binary | std::ios::in);
    if (!in)
    {
        std::cerr << "Does not exist." << std::endl;
        return 1;
    }

    // Load sbi.numAlts from file.
    in.read((char*)&sbi.numAlts,sizeof(sbi.numAlts));
    sbi.numAltsBase = std::log(sbi.numAlts) / std::log(2);

    // Load the barcode table of the index from file.
    sbi.barcode_table.load(in);
    in.close();

    sbi.rank_support_barcode_table = sdsl::rank_support_v<>(&sbi.barcode_table);

    // Load the match table of the index from file.
    CharString matchFilename = filename;
    matchFilename += ".match";
    bool ret = sdsl::load_from_file(sbi.match_table, toCString(matchFilename));
    if (ret == false)
    {
        std::cerr << "Match table does not exist." << std::endl;
        return 2;
    }

    sbi.rank_support_match_table = sdsl::rank_support_v<>(&sbi.match_table);

    // Load the substitution table of the index from file.
    CharString substFilename = filename;
    substFilename += ".subst";
    ret = sdsl::load_from_file(sbi.substitution_table, toCString(substFilename));
    if (ret == false)
    {
        std::cerr << "Substitution table does not exist." << std::endl;
        return 3;
    }

    printDone();
    return 0;
}

void write(seqan::CharString & filename, BarcodeIndex & sbi)
{
    writeBarcodeTable(filename, sbi);
    writeMatchTable(filename, sbi);
    writeSubstitutionTable(filename, sbi);
}

void writeBarcodeTable(seqan::CharString & filename, BarcodeIndex & sbi)
{
    printStatus("Writing barcode table of barcode index");

    CharString bcFilename = filename;
    bcFilename += ".bc";

    // Open the barcode table file.
    std::ofstream out;
    out.open(toCString(bcFilename), std::ios::binary | std::ios::trunc | std::ios::out);
    if (!out)
    {
        std::cerr << "Cannot open file '" << toCString(bcFilename) << "'." << std::endl;
        return;
    }

    // Write sbi.numAlts to file.
    out.write((char*)&sbi.numAlts,sizeof(sbi.numAlts));

    // Write the barcode table of the index to file.
    sbi.barcode_table.serialize(out);
    out.close();

    printDone();
}

void writeMatchTable(seqan::CharString & filename, BarcodeIndex & sbi)
{
    printStatus("Writing match table of barcode index");

    CharString matchFilename = filename;
    matchFilename += ".match";

    // Write the match table of the index to a file.
    sdsl::store_to_file(sbi.match_table, toCString(matchFilename));

    printDone();
}

void writeSubstitutionTable(seqan::CharString & filename, BarcodeIndex & sbi)
{
    printStatus("Writing substitution table of barcode index");

    CharString substFilename = filename;
    substFilename += ".subst";

    // Write the substitution table of the index to a file.
    sdsl::store_to_file(sbi.substitution_table, toCString(substFilename));

    printDone();
}

BarcodeStatus retrieve(std::vector<seqan::DnaString> & bx, BarcodeIndex & sbi, seqan::DnaString & rx, seqan::CharString & qx)
{
    uint64_t h = hash(rx);
    BarcodeStatus s = get_status(sbi, h);
    std::vector<std::pair<DnaString, unsigned> > bxx;

    switch (s)
    {
        case BarcodeStatus::UNRECOGNIZED:
            break;
        case BarcodeStatus::INVALID:
            break;
        case BarcodeStatus::MATCH:
        {
            bx.push_back(rx);
            break;
        }
        case BarcodeStatus::ONE_ERROR:
        {
            unsigned offset = 0;
            unsigned i;
            while(offset != sbi.numAlts && get_substitution(i, sbi, h, offset))
            {
                uint64_t h_corrected = get_corrected_barcode(sbi, h, i);
                SEQAN_ASSERT_NEQ(h_corrected, h);
                bxx.push_back(std::pair<DnaString, unsigned>(unhash(h_corrected, sbi.bcLength), qx[length(qx)-1 - i]));
                ++offset;
            }

            std::sort(bxx.begin(), bxx.end(), [](auto & left, auto & right) {
                return left.second < right.second;
            });

            for (unsigned i = 0; i < bxx.size(); ++i)
                bx.push_back(bxx[i].first);

            break;
        }
    }

    return s;
}

BarcodeStatus retrieve(std::vector<seqan::DnaString> & bx, BarcodeIndex & sbi, seqan::Dna5String & rx, seqan::CharString & qx)
{
    // Find N positions in the barcode.
    typename Iterator<Dna5String, Rooted>::Type it = begin(rx);
    typename Iterator<Dna5String, Rooted>::Type itEnd = end(rx);
    String<unsigned> posN;
    for (; it != itEnd; ++it)
    {
        if (*it == 'N')
            appendValue(posN, position(it));
    }

    // Call retrieve() for each possible nucleotide at N positions in the barcode.
    DnaString rxx = rx;
    switch(length(posN))
    {
        case 0:
            return retrieve(bx, sbi, rxx, qx);
        case 1:
        {
            // Determine whether the barcode status will be ONE_ERROR or UNRECOGNIZED.
            BarcodeStatus ret = BarcodeStatus::UNRECOGNIZED;
            for (unsigned i = 0; i < ValueSize<Dna>::VALUE; ++i)
            {
                rxx[posN[0]] = Dna(i);
                uint64_t h = hash(rxx);
                BarcodeStatus s = get_status(sbi, h);
                if (s == BarcodeStatus::MATCH)
                {
                    // Add N substituted MATCH barcode as ONE_ERROR barcode.
                    ret = BarcodeStatus::ONE_ERROR;
                    bx.push_back(rxx);
                }
            }
            return ret;
        }
        default: // more than 2 Ns in the barcode
            return BarcodeStatus::UNRECOGNIZED;
    }
    return BarcodeStatus::INVALID;
}
