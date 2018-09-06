#ifndef BARCODE_INDEX_H_
#define BARCODE_INDEX_H_

#include <sdsl/bit_vectors.hpp>
#include <seqan/sequence.h>

// -----------------------------------------------------------------------------
// Succinct barcode index
// -----------------------------------------------------------------------------


struct BarcodeIndex
{
    unsigned numAlts;
    unsigned numAltsBase;
    unsigned bcLength;
    sdsl::bit_vector barcode_table;
    sdsl::rank_support_v<> rank_support_barcode_table;
    sdsl::bit_vector match_table;
    sdsl::rank_support_v<> rank_support_match_table;
    sdsl::int_vector<> substitution_table;
    
    BarcodeIndex() {}    
    BarcodeIndex(seqan::CharString & filename);
};

enum class BarcodeStatus {
    UNRECOGNIZED,
    INVALID,
    MATCH,
    ONE_ERROR
};

void buildIndex(BarcodeIndex & sbi, seqan::CharString & filename, unsigned alternatives);
void buildSubstitutionTable(BarcodeIndex & sbi, seqan::CharString & filename);
int load(BarcodeIndex & sbi, seqan::CharString & filename);
void write(seqan::CharString & filename, BarcodeIndex & sbi);
void writeBarcodeTable(seqan::CharString & filename, BarcodeIndex & sbi);
void writeMatchTable(seqan::CharString & filename, BarcodeIndex & sbi);
void writeSubstitutionTable(seqan::CharString & filename, BarcodeIndex & sbi);
BarcodeStatus retrieve(std::vector<seqan::DnaString> & bx, BarcodeIndex & sbi, seqan::DnaString & rx, seqan::CharString & qx);
BarcodeStatus retrieve(std::vector<seqan::DnaString> & bx, BarcodeIndex & sbi, seqan::Dna5String & rx, seqan::CharString & qx);

#endif // BARCODE_INDEX_H_
