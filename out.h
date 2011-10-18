#ifndef out_h
#define out_h

#include "seq.h"
#include "hsp.h"
#include "orf.h"
#include "path.h"
#include <float.h>

struct AlnData
	{
/***
SA.Seq and SB.Seq align.
Reverse strand stuff for nucleotides is handled like this:
	SA.RevComp must be false.
	If SB.RevComp is true, then SA.Seq is r.c.'d relative to the sequence in
	the input file (query or db). If so, coordinates in HSP refer to SB.Seq
	so are also r.c.'d relative to the original sequence.
***/
	SeqData SA;
	SeqData SB;
	HSPData HSP;
	const char *Path;
	char IdDesc[256];

	float FractId;
	float RawScore;
	float BitScore;
	float Evalue;

	void LogMe() const
		{
		Log("AD: ");
		HSP.LogMe();
		Log(" %s,%s\n", SA.Label, SB.Label);
		}
	};

bool OnDerepHit(const SeqData &SA, const SeqData &SB);

bool OnLocalUngappedHit(const SeqData &SA, const SeqData &SB,
  const HSPData &HSP, float &Evalue, float &FractId);

bool OnLocalGappedHit(const SeqData &SA, const SeqData &SB,
  const HSPData &HSP, const PathData &PD, float &Evalue, float &FractId);

bool OnGlobalHit(const SeqData &SA, const SeqData &SB, const PathData &PD,
  float &FractId);

void OnReject(const SeqData &SA, const SeqData &SB, double FractId,
  const char *Path);

void OnNotMatched(const char *Label, unsigned L);
void OnNewCluster(unsigned ClusterIndex, const char *Label, unsigned L);
void OnNewLibCluster(unsigned ClusterIndex, const char *Label, unsigned L);
void OnLibCluster(unsigned ClusterIndex, unsigned Size, double AvgId,
  const char *Label);
void OnNewCluster(unsigned ClusterIndex, unsigned Size, double AvgId,
  const char *Label);
void OnChainCov(const SeqData &NucleoSD, const SeqData &TargetSD,
  float Score, float ChainCov);

void SetUserFieldIndexes(const string &s);

void BlastOut(FILE *f, const AlnData &AD);
void Blast6Out(FILE *f, const AlnData &AD);
void FastaPairOut(FILE *f, const AlnData &AD);
void UserOut(FILE *f, const AlnData &AD);

void BlastOutORF(FILE *f, const AlnData &AD);

void OpenOutputFiles();
void CloseOutputFiles();
void SetLibSeedCount(unsigned DBSeqCount);
const char *UserFieldIndexToStr(unsigned i);

extern float **g_SubstMx;

static char g_IdChar = '|';
static char g_DiffChar = ' ';

static inline char GetSymN(byte Letter1, byte Letter2)
	{
	Letter1 = toupper(Letter1);
	Letter2 = toupper(Letter2);
	if (Letter1 == Letter2)
		return g_IdChar;
	return g_DiffChar;
	}

static inline char GetSymA(byte Letter1, byte Letter2)
	{
	Letter1 = toupper(Letter1);
	Letter2 = toupper(Letter2);
	if (Letter1 == Letter2)
		return '|';

	float Score = g_SubstMx[Letter1][Letter2];
	if (Score >= 2.0f)
		return ':';
	if (Score > 0.0f)
		return '.';
	return ' ';
	}

static inline char GetSym(byte Letter1, byte Letter2, bool Nucleo)
	{
	if (Nucleo)
		return GetSymN(Letter1, Letter2);
	else
		return GetSymA(Letter1, Letter2);
	}

static unsigned GetNDig(unsigned n)
	{
	if (n < 10)
		return 1;
	if (n < 100)
		return 2;
	if (n < 1000)
		return 3;
	if (n < 10000)
		return 4;
	if (n < 100000)
		return 5;
	if (n < 1000000)
		return 6;
	return 10;
	}

extern unsigned *g_UserFieldIndexes;
extern unsigned g_UserFieldCount;

#endif // out_h
