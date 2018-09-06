#ifndef UTILS_H_
#define UTILS_H_

#include <sstream>
#include <seqan/sequence.h>

uint64_t hash(seqan::DnaString & barcode);
seqan::DnaString unhash(uint64_t h, unsigned bcLength);

std::string currentTime();

void printStatus(const char * message);
void printStatus(std::ostringstream & message);

void printDone();

void printInfo(const char * message);
void printInfo(std::ostringstream & message);

void printWarning(const char * message);
void printWarning(std::ostringstream & message);

#endif  // UTILS_H_