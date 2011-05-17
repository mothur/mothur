//uchime by Robert C. Edgar http://drive5.com/uchime This code is donated to the public domain.

#ifndef alpha_h
#define alpha_h

#include <limits.h>
#include <string>

using namespace std;

const unsigned INVALID_LETTER = 0;
const unsigned char INVALID_CHAR = '?';

extern unsigned g_CharToLetterAmino[];
extern unsigned g_CharToLetterAminoStop[];
extern unsigned char g_LetterToCharAmino[];
extern unsigned g_CharToLetterNucleo[];
extern unsigned char g_LetterToCharNucleo[];
extern unsigned g_CodonWordToAminoLetter[];
extern char g_CodonWordToAminoChar[];
extern unsigned char g_CharToCompChar[];
extern unsigned g_CharToCompLetter[];
extern bool g_IsAminoChar[];
extern bool g_IsNucleoChar[];
extern bool g_IsACGTU[];
extern float g_AminoFreqs[];

extern unsigned g_CharToLetterRed[];
extern unsigned char g_LetterToCharRed[];
extern unsigned g_RedAlphaSize;

void LogRedAlphaRed();
void ReadRedAlphaFromFile(const string &FileName);
unsigned char GetAminoCharFrom3NucChars(unsigned char c1, unsigned char c2,
  unsigned char c3);

static inline bool AminoLetterIsStartCodon(unsigned char Letter)
	{
	return Letter == 10;
	}

static inline bool AminoLetterIsStopCodon(unsigned char Letter)
	{
	return Letter == 20;
	}

const char *WordToStr(unsigned Word, unsigned WordLength, bool Nucleo);
const char *WordToStrNucleo(unsigned Word, unsigned WordLength);
const char *WordToStrAmino(unsigned Word, unsigned WordLength);
const char *WordToStrAmino2(unsigned Word, unsigned WordLength, char *Str);

#endif // alpha_h
