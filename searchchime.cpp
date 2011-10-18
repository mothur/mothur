#include "myutils.h"
#include "ultra.h"
#include "chime.h"
#include "uc.h"
#include "dp.h"
#include <set>
#include <algorithm>

#define TRACE	0

extern FILE *g_fUChime;

void GetCandidateParents(Ultra &U, const SeqData &QSD, float AbQ,
  vector<unsigned> &Parents);

void AlignChime(const SeqData &QSD, const SeqData &ASD, const SeqData &BSD,
  const string &PathQA, const string &PathQB, ChimeHit2 &Hit);

double GetFractIdGivenPath(const byte *A, const byte *B, const char *Path, bool Nucleo);

static void GetSmoothedIdVec(const SeqData &QSD, const SeqData &PSD, const string &Path,
  vector<unsigned> &IdVec, unsigned d)
	{
	IdVec.clear();
	const unsigned ColCount = SIZE(Path);

	const byte *Q = QSD.Seq;
	const byte *P = PSD.Seq;

	const unsigned QL = QSD.L;
	const unsigned PL = PSD.L;

	if (QL <= d)
		{
		IdVec.resize(QSD.L, 0);
		return;
		}

	unsigned QPos = 0;
	unsigned PPos = 0;

	vector<bool> SameVec;
	SameVec.reserve(QL);
	for (unsigned Col = 0; Col < ColCount; ++Col)
		{
		char c = Path[Col];

		bool Same = false;
		if (c == 'M')
			{
			byte q = Q[QPos];
			byte p = P[PPos];
			Same = (toupper(q) == toupper(p));
			}

		if (c == 'M' || c == 'D')
			{
			++QPos;
			SameVec.push_back(Same);
			}

		if (c == 'M' || c == 'I')
			++PPos;
		}

	asserta(SIZE(SameVec) == QL);

	unsigned n = 0;
	for (unsigned QPos = 0; QPos < d; ++QPos)
		{
		if (SameVec[QPos])
			++n;
		IdVec.push_back(n);
		}

	for (unsigned QPos = d; QPos < QL; ++QPos)
		{
		if (SameVec[QPos])
			++n;
		IdVec.push_back(n);
		if (SameVec[QPos-d])
			--n;
		}
	asserta(SIZE(IdVec) == QL);

#if	TRACE
	{
	Log("\n");
	Log("GetSmoothedIdVec\n");
	unsigned QPos = 0;
	unsigned PPos = 0;
	Log("Q P  Same       Id\n");
	Log("- -  ----  -------\n");
	for (unsigned Col = 0; Col < ColCount; ++Col)
		{
		char c = Path[Col];

		bool Same = false;
		if (c == 'M')
			{
			byte q = Q[QPos];
			byte p = P[PPos];
			Same = (toupper(q) == toupper(p));
			Log("%c %c  %4c  %7d\n", q, p, tof(Same), IdVec[QPos]);
			}

		if (c == 'M' || c == 'D')
			++QPos;
		if (c == 'M' || c == 'I')
			++PPos;
		}
	}
#endif
	}

bool SearchChime(Ultra &U, const SeqData &QSD, float QAb, 
  const AlnParams &AP, const AlnHeuristics &AH, HSPFinder &HF,
  float MinFractId, ChimeHit2 &Hit)
	{
	Hit.Clear();
	Hit.QLabel = QSD.Label;

	if (opt_verbose)
		{
		Log("\n");
		Log("SearchChime()\n");
		Log("Query>%s\n", QSD.Label);
		}

	vector<unsigned> Parents;
	GetCandidateParents(U, QSD, QAb, Parents);

	unsigned ParentCount = SIZE(Parents);
	if (ParentCount <= 1)
		{
		if (opt_verbose)
			Log("%u candidate parents, done.\n", ParentCount);
		return false;
		}

	if (opt_fastalign)
		HF.SetA(QSD);
	HSPFinder *ptrHF = (opt_fastalign ? &HF : 0);

	unsigned ChunkLength;
	vector<unsigned> ChunkLos;
	GetChunkInfo(QSD.L, ChunkLength, ChunkLos);
	const unsigned ChunkCount = SIZE(ChunkLos);

	vector<unsigned> ChunkIndexToBestId(ChunkCount, 0);
	vector<unsigned> ChunkIndexToBestParentIndex(ChunkCount, UINT_MAX);

	vector<SeqData> PSDs;
	vector<string> Paths;
	double TopPctId = 0.0;
	unsigned TopParentIndex = UINT_MAX;
	unsigned QL = QSD.L;
	vector<unsigned> MaxIdVec(QL, 0);
	for (unsigned ParentIndex = 0; ParentIndex < ParentCount; ++ParentIndex)
		{
		unsigned ParentSeqIndex = Parents[ParentIndex];

		SeqData PSD;
		//PSD.Label = U.GetSeedLabel(ParentSeqIndex);
		//PSD.Seq = U.GetSeedSeq(ParentSeqIndex);
		//PSD.L = U.GetSeedLength(ParentSeqIndex);
		//PSD.Index = ParentSeqIndex;
		U.GetSeqData(ParentSeqIndex, PSD);
		PSDs.push_back(PSD);

		if (opt_fastalign)
			HF.SetB(PSD);

		PathData PD;

		float HSPId;
		bool Found = GlobalAlign(QSD, PSD, AP, AH, *ptrHF, MinFractId, HSPId, PD);
		if (!Found)
			{
			Paths.push_back("");				
			continue;
			}

		double PctId = 100.0*GetFractIdGivenPath(QSD.Seq, PSD.Seq, PD.Start, true);
		if (opt_selfid && PctId == 100.0)
			{
			Paths.push_back("");				
			continue;
			}

		if (PctId > TopPctId)
			{
			TopParentIndex = ParentIndex;
			TopPctId = PctId;
			if (TopPctId >= 100.0 - opt_mindiv)
				{
				if (opt_verbose)
					{
					Log("  %.1f%%  >%s\n", TopPctId, PSD.Label);
					Log("  Top hit exceeds ctl threshold, done.\n");
					return false;
					}
				}
			}

		string Path = PD.Start;
		Paths.push_back(Path);

		vector<unsigned> IdVec;
		GetSmoothedIdVec(QSD, PSD, Path, IdVec, opt_idsmoothwindow);

		for (unsigned QPos = 0; QPos < QL; ++QPos)
			if (IdVec[QPos] > MaxIdVec[QPos])
				MaxIdVec[QPos] = IdVec[QPos];
		}

	vector<unsigned> BestParents;
	for (unsigned k = 0; k < opt_maxp; ++k)
		{
		unsigned BestParent = UINT_MAX;
		unsigned BestCov = 0;
		for (unsigned ParentIndex = 0; ParentIndex < ParentCount; ++ParentIndex)
			{
			const SeqData &PSD = PSDs[ParentIndex];
			const string &Path = Paths[ParentIndex];
			if (Path == "")
				continue;

			vector<unsigned> IdVec;
			GetSmoothedIdVec(QSD, PSD, Path, IdVec, opt_idsmoothwindow);

			unsigned Cov = 0;
			for (unsigned QPos = 0; QPos < QL; ++QPos)
				if (IdVec[QPos] == MaxIdVec[QPos])
					++Cov;

			if (Cov > BestCov)
				{
				BestParent = ParentIndex;
				BestCov = Cov;
				}
			}

		if (BestParent == UINT_MAX)
			break;

		BestParents.push_back(BestParent);
		vector<unsigned> IdVec;

		const SeqData &PSD = PSDs[BestParent];
		const string &Path = Paths[BestParent];
		GetSmoothedIdVec(QSD, PSD, Path, IdVec, opt_idsmoothwindow);
		for (unsigned QPos = 0; QPos < QL; ++QPos)
			if (IdVec[QPos] == MaxIdVec[QPos])
				MaxIdVec[QPos] = UINT_MAX;
		}

	unsigned BestParentCount = SIZE(BestParents);

	if (opt_verbose)
		{
		Log("%u/%u best parents\n", BestParentCount, ParentCount);
		for (unsigned k = 0; k < BestParentCount; ++k)
			{
			unsigned i = BestParents[k];
			Log(" %s\n", PSDs[i].Label);
			}
		}

	bool Found = false;
	for (unsigned k1 = 0; k1 < BestParentCount; ++k1)
		{
		unsigned i1 = BestParents[k1];
		asserta(i1 < ParentCount);

		const SeqData &PSD1 = PSDs[i1];
		const string &Path1 = Paths[i1];

		for (unsigned k2 = k1 + 1; k2 < BestParentCount; ++k2)
			{
			unsigned i2 = BestParents[k2];
			asserta(i2 < ParentCount);
			asserta(i2 != i1);

			const SeqData &PSD2 = PSDs[i2];
			const string &Path2 = Paths[i2];

			ChimeHit2 Hit2;
			AlignChime(QSD, PSD1, PSD2, Path1, Path2, Hit2);
			Hit2.PctIdQT = TopPctId;

			if (Hit2.Accept())
				Found = true;

			if (Hit2.Score > Hit.Score)
				Hit = Hit2;

			if (opt_verbose)
				Hit2.LogMe();
			}
		}

	return Found;
	}
