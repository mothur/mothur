//uchime by Robert C. Edgar http://drive5.com/uchime This code is donated to the public domain.

#include "myutils.h"
#include "alpha.h"
#include "timing.h"

bool isgap(byte c)
	{
	return c == '-' || c == '.';
	}

const char *WordToStrAmino(unsigned Word, unsigned WordLength)
	{
	static char Str[32];
	for (unsigned i = 0; i < WordLength; ++i)
		{
		unsigned Letter = Word%20;
		Str[WordLength-i-1] = g_LetterToCharAmino[Letter];
		Word /= 20;
		}
	Str[WordLength] = 0;
	return Str;
	}

const char *WordToStrAmino2(unsigned Word, unsigned WordLength, char *Str)
	{
	for (unsigned i = 0; i < WordLength; ++i)
		{
		unsigned Letter = Word%20;
		Str[WordLength-i-1] = g_LetterToCharAmino[Letter];
		Word /= 20;
		}
	Str[WordLength] = 0;
	return Str;
	}

const char *WordToStrNucleo(unsigned Word, unsigned WordLength)
	{
	static char Str[32];
	for (unsigned i = 0; i < WordLength; ++i)
		{
		unsigned Letter = Word%4;
		Str[WordLength-i-1] = g_LetterToCharNucleo[Letter];
		Word /= 4;
		}
	Str[WordLength] = 0;
	return Str;
	}

const char *WordToStr(unsigned Word, unsigned WordLength, bool Nucleo)
	{
	return (Nucleo ? WordToStrNucleo : WordToStrAmino)(Word, WordLength);
	}

byte *RevCompAlloc(const byte *Seq, unsigned L)
	{
	byte *RCSeq = MYALLOC(byte, L, Alpha);

	for (unsigned i = 0; i < L; ++i)
		RCSeq[L-i-1] = g_CharToCompChar[Seq[i]];

	return RCSeq;
	}

void RevCompInPlace(byte *Seq, unsigned L)
	{
	unsigned L1 = L - 1;
	unsigned L2 = L/2;
	for (unsigned i = 0; i < L2; ++i)
		{
		unsigned j = L1 - i;
		unsigned ci = Seq[i];
		unsigned cj = Seq[j];

		unsigned ri = g_CharToCompChar[ci];
		unsigned rj = g_CharToCompChar[cj];

		Seq[i] = rj;
		Seq[j] = ri;
		}

	if (L%2 == 1)
		Seq[L2] = g_CharToCompChar[Seq[L2]];
	}

void RevComp(const byte *Seq, unsigned L, byte *RCSeq)
	{
	for (unsigned i = 0; i < L; ++i)
		RCSeq[L-i-1] = g_CharToCompChar[Seq[i]];
	}

unsigned char GetAminoCharFrom3NucChars(unsigned char c1, unsigned char c2,
  unsigned char c3)
	{
	unsigned Letter1 = g_CharToLetterNucleo[c1];
	unsigned Letter2 = g_CharToLetterNucleo[c2];
	unsigned Letter3 = g_CharToLetterNucleo[c3];
	unsigned Word = Letter1*(4*4) + Letter2*4 + Letter3;

	unsigned Letter = g_CodonWordToAminoLetter[Word];
	return g_LetterToCharAmino[Letter];
	}
