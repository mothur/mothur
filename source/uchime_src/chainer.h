#ifndef chainer_h
#define chainer_h

#include "hsp.h"
#include "seq.h"
#include <list>

const float BAD_SCORE = -9e9f;

struct TargetHit
	{
	unsigned TargetIndex;
	unsigned TargetLo;
	unsigned TargetHi;
	int QueryFrame;
	float RawScore; // SOMETIMES USED FOR BIT SCORE!!!
//	unsigned TargetLength;

	void LogMe() const
		{
		Log("lo %u, hi %u, frame %d, score %.1f\n",
		  TargetLo, TargetHi, QueryFrame, RawScore);
		}
	};

struct ChainData
	{
	unsigned LastHSPIndex;
	unsigned Ahi;
	unsigned Bhi;
	float Score;
	};

class Chainer
	{
public:
	HSPData **m_HSPs; // memory owned elsewhere
	unsigned m_HSPCount;
	unsigned m_MaxHSPCount;

	BPData *m_BPs;

	unsigned *m_PrevHSPIndexes;		// Predecessor in chain
	float *m_HSPIndexToChainScore;

	list<unsigned> m_Chains;		// Live HSP indexes

public:
	Chainer();
	~Chainer();
	void Reset();
	void Clear(bool ctor = false);
	float Chain(HSPData **HSPs, unsigned HSPCount, HSPData **OptChain,
	  unsigned &OptChainLength);
	bool ResolveOverlaps(const SeqData &SA, const SeqData &SB, double MinScore,
	  const float * const *SubstMx, HSPData **InHSPs, unsigned InHSPCount,
	  HSPData **OutHSPs, unsigned &OutHSPCount);
	void ResolveOverlap(HSPData &HSP1, HSPData &HSP2);

	float ChainBrute(HSPData **HSPs, unsigned HSPCount, HSPData **OptChain,
	  unsigned &OptChainLength);
	void LogMe() const;
	void LogHSPs(HSPData **HSPs, unsigned HSPCount) const;
	void LogBPs() const;

	static bool IsValidChain(HSPData **HSPs, unsigned HSPCount);
	static void AssertValidChain(HSPData **HSPs, unsigned HSPCount);
	static void LogChain(HSPData **HSPs, unsigned HSPCount);
	static void LogChain2(HSPData **HSPs, unsigned HSPCount);
	static float GetChainScore(HSPData **HSPs, unsigned HSPCount);

private:
	void AllocHSPCount(unsigned MaxHSPCount);
	void SetBPs();
	void SortBPs();
	unsigned FindBestChainLT(unsigned Ahi, unsigned Bhi);
	};

#endif // chainer_h
