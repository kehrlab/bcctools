#ifndef UTILS_H_
#define UTILS_H_

#include <sstream>
#include <seqan/sequence.h>

uint64_t hash(seqan::DnaString & barcode);
seqan::DnaString unhash(uint64_t h, unsigned bcLength);

bool union_by_index(std::vector<unsigned> & uf, unsigned a, unsigned b);
unsigned find(std::vector<unsigned> & uf, unsigned a);

std::vector<seqan::DnaString> parseBarcodeList(const char * cBarcode);

std::string currentTime();

void printStatus(const char * message);
void printStatus(std::ostringstream & message);

void printDone();

void printInfo(const char * message);
void printInfo(std::ostringstream & message);

void printWarning(const char * message);
void printWarning(std::ostringstream & message);

#endif  // UTILS_H_