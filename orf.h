#ifndef orf_h
#define orf_h

#include "alpha.h"

struct ORFData
	{
	const byte *NucSeq;
	const byte *AminoSeq;
	int Frame;
	unsigned NucL;
	unsigned AminoL;
	unsigned NucLo;
	unsigned NucHi;
	ORFData *Next;

	unsigned GetNucPosFirstBase() const;
	unsigned GetAAPos(unsigned NucPos) const;
	unsigned GetCodex(unsigned NucPos) const;
	unsigned GetNucLo(unsigned AALo, unsigned AAHi) const;
	unsigned GetNucHi(unsigned AALo, unsigned AAHi) const;
	unsigned GetAALo(unsigned NucLo, unsigned NucHi) const;
	unsigned GetAAHi(unsigned NucLo, unsigned NucHi) const;
	unsigned GetNucPosFirstBaseInCodon(unsigned AAPos) const;
	unsigned GetNucPosLastBaseInCodon(unsigned AAPos) const;
	unsigned RoundToCodonLo(unsigned NucPos) const;
	unsigned RoundToCodonHi(unsigned NucPos) const;
	void LogMe() const;
	void LogMe2() const;
	};

const byte ORFEND = '.';

void GetORFs(const byte *NucSeq, unsigned NucL, vector<ORFData> &ORFs,
  unsigned ORFStyle, int FindFrame, int Sign);

#endif // orf_h
