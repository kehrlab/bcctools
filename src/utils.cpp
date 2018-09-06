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
    std::cerr << "[bctools " << currentTime() << "] " << message << "... " << std::flush;
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
    std::cerr << "[bctools " << currentTime() << "] " << message << std::endl;
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
    std::cerr << "[bctools WARNING] " << message << std::endl;
}

void printWarning(std::ostringstream & message)
{
    std::string msg = message.str();
    printWarning(toCString(msg));
    message.str("");
}