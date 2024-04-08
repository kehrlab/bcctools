#include <sstream>
#include <ctime>

#include <seqan/sequence.h>

using namespace seqan;

// ---------------------------------------------------------------------------------------
// Function hash()
// ---------------------------------------------------------------------------------------

uint64_t hash(DnaString & barcode)
{
    uint64_t h = 0;
    Iterator<DnaString>::Type it = begin(barcode);
    Iterator<DnaString>::Type itEnd = end(barcode);
    for (; it != itEnd; ++it)
        h = (h << 2) | ordValue((Dna)*it); // Add i^th position to hash value:
                             // Multiply by 4, add ordValue of next base pair.
    return h;
}

// ---------------------------------------------------------------------------------------
// Function unhash()
// ---------------------------------------------------------------------------------------

DnaString unhash(uint64_t h, unsigned bcLength)
{
    DnaString barcode;
    resize(barcode, bcLength);
    for (unsigned i = bcLength; i > 0; )
    {
        barcode[--i] = (Dna)(h % ValueSize<Dna>::VALUE);
        h /= ValueSize<Dna>::VALUE;
    }
    return barcode;
}

// ---------------------------------------------------------------------------------------
// Functions union() and find()
// ---------------------------------------------------------------------------------------

bool union_by_index(std::vector<unsigned> & uf, unsigned a, unsigned b)
{
    if (a < uf.size() && b < uf.size())
    {
        uf[b] = uf[a];
        return true;
    }
    return false;
}

unsigned find(std::vector<unsigned> & uf, unsigned a)
{
    while (uf[a] != a)
    {
        uf[a] = uf[uf[a]];
        a = uf[a];
    }
    return a;
}

// ---------------------------------------------------------------------------------------
// Function parseBarcodeList()
// ---------------------------------------------------------------------------------------

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

// ---------------------------------------------------------------------------------------
// Function currentTime()
// ---------------------------------------------------------------------------------------

std::string currentTime()
{
    // Get the current date and time.
    time_t now = time(0);
    struct tm tstruct;
    tstruct = *localtime(&now);

    // Format the time to 'YYYY-MM-DD HH:MM:SS'.
    char timeStamp[80];
    strftime(timeStamp, sizeof(timeStamp), "%Y-%m-%d %X", &tstruct);

    return std::string(timeStamp);
}

// ---------------------------------------------------------------------------------------
// Function printStatus()
// ---------------------------------------------------------------------------------------

void printStatus(const char * message)
{
    // Print time and message *without* newline.
    std::cerr << "[bcctools " << currentTime() << "] " << message << "... " << std::flush;
}

void printStatus(std::ostringstream & message)
{
    std::string msg = message.str();
    printStatus(toCString(msg));
    message.str("");
}

// ---------------------------------------------------------------------------------------
// Function printDone()
// ---------------------------------------------------------------------------------------

void printDone()
{
    // Print 'Done.' and newline.
    std::cerr << "Done." << std::endl;
}

// ---------------------------------------------------------------------------------------
// Function printInfo()
// ---------------------------------------------------------------------------------------

void printInfo(const char * message)
{
    // Print time, message and newline.
    std::cerr << "[bcctools " << currentTime() << "] " << message << std::endl;
}

void printInfo(std::ostringstream & message)
{
    std::string msg = message.str();
    printInfo(toCString(msg));
    message.str("");
}

// ---------------------------------------------------------------------------------------
// Function printWarning()
// ---------------------------------------------------------------------------------------

void printWarning(const char * message)
{
    // Print warning message and newline.
    std::cerr << "[bcctools WARNING] " << message << std::endl;
}

void printWarning(std::ostringstream & message)
{
    std::string msg = message.str();
    printWarning(toCString(msg));
    message.str("");
}