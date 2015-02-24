#include "myutils.h"
#include "seq.h"
#include "chime.h"

#define	TRACE	0

/***
Let:
	S[i] =	Score of col i: 0=no SNP, +1 = Y, -3 = N or A.

	V[k] =	Best segment score from j, j+1 .. k for all possible j
			max(j) Sum i=j..k S[i]

Recursion relation:
	V[k] =	S[k] + max (V[k-1], 0)
***/

void AlignChimeGlobal3(const string &Q3, const string &A3, const string &B3,
  const string &QLabel, const string &ALabel, const string &BLabel,
  ChimeHit2 &Hit);

void Make3Way(const SeqData &SDQ, const SeqData &SDA, const SeqData &SDB,
  const string &PathQA, const string &PathQB,
  string &Q3, string &A3, string &B3);

double GetScore2(double Y, double N, double A);

void AlignChimeLocal3(const string &Q3, const string &A3, const string &B3,
  const string &QLabel, const string &ALabel, const string &BLabel,
  ChimeHit2 &Hit)
	{
	Hit.Clear();

	const byte *Q3Seq = (const byte *) Q3.c_str();
	const byte *A3Seq = (const byte *) A3.c_str();
	const byte *B3Seq = (const byte *) B3.c_str();

	const unsigned ColCount = SIZE(Q3);
	asserta(SIZE(A3) == ColCount && SIZE(B3) == ColCount);

	vector<float> ColScoresA(ColCount, 0.0f);
	vector<float> ColScoresB(ColCount, 0.0f);

	float ScoreN = -(float) opt_xn;
	unsigned QL = 0;
	for (unsigned Col = 0; Col < ColCount; ++Col)
		{
		char q = Q3Seq[Col];
		char a = A3Seq[Col];
		char b = B3Seq[Col];

		if (!isgap(q))
			++QL;

		if (q == a && q == b && a == b)
			continue;

		if (isgap(q) || isgap(a) || isgap(b))
			continue;

		if (Col > 0 && (isgap(Q3Seq[Col-1]) || isgap(A3Seq[Col-1]) || isgap(B3Seq[Col-1])))
			continue;

		if (Col + 1 < ColCount && (isgap(Q3Seq[Col+1]) || isgap(A3Seq[Col+1]) || isgap(B3Seq[Col+1])))
			continue;

		if (q == a && q != b)
			ColScoresA[Col] = 1;
		else
			ColScoresA[Col] = ScoreN;

		if (q == b && q != a)
			ColScoresB[Col] = 1;
		else
			ColScoresB[Col] = ScoreN;
		}

	vector<float> LVA(ColCount, 0.0f);
	vector<float> LVB(ColCount, 0.0f);

	LVA[0] = ColScoresA[0];
	LVB[0] = ColScoresB[0];
	for (unsigned Col = 1; Col < ColCount; ++Col)
		{
		LVA[Col] = max(LVA[Col-1], 0.0f) + ColScoresA[Col];
		LVB[Col] = max(LVB[Col-1], 0.0f) + ColScoresB[Col];
		}

	vector<float> RVA(ColCount, 0.0f);
	vector<float> RVB(ColCount, 0.0f);

	RVA[ColCount-1] = ColScoresA[ColCount-1];
	RVB[ColCount-1] = ColScoresB[ColCount-1];
	for (int Col = ColCount-2; Col >= 0; --Col)
		{
		RVA[Col] = max(RVA[Col+1], 0.0f) + ColScoresA[Col];
		RVB[Col] = max(RVB[Col+1], 0.0f) + ColScoresB[Col];
		}

	bool FirstA = true;
	float MaxSum = 0.0;
	unsigned ColX = UINT_MAX;
	for (unsigned Col = 1; Col < ColCount-1; ++Col)
		{
		float Sum = LVA[Col] + RVB[Col+1];
		if (Sum > MaxSum)
			{
			FirstA = true;
			MaxSum = Sum;
			ColX = Col;
			}
		}

	for (unsigned Col = 1; Col < ColCount-1; ++Col)
		{
		float Sum = LVB[Col] + RVA[Col+1];
		if (Sum > MaxSum)
			{
			FirstA = false;
			MaxSum = Sum;
			ColX = Col;
			}
		}
	if (ColX == UINT_MAX)
		return;

	unsigned ColLo = UINT_MAX;
	unsigned ColHi = UINT_MAX;
	if (FirstA)
		{
		float Sum = 0.0f;
		for (int Col = ColX; Col >= 0; --Col)
			{
			Sum += ColScoresA[Col];
			if (Sum >= LVA[ColX])
				{
				ColLo = Col;
				break;
				}
			}
		asserta(Sum >= LVA[ColX]);
		Sum = 0.0f;
		for (unsigned Col = ColX+1; Col < ColCount; ++Col)
			{
			Sum += ColScoresB[Col];
			if (Sum >= RVB[ColX])
				{
				ColHi = Col;
				break;
				}
			}
		asserta(Sum >= RVB[ColX]);
		}
	else
		{
		float Sum = 0.0f;
		for (int Col = ColX; Col >= 0; --Col)
			{
			Sum += ColScoresB[Col];
			if (Sum >= LVB[ColX])
				{
				ColLo = Col;
				break;
				}
			}
		asserta(Sum >= LVB[ColX]);
		Sum = 0.0f;
		for (unsigned Col = ColX+1; Col < ColCount; ++Col)
			{
			Sum += ColScoresA[Col];
			if (Sum >= RVA[ColX])
				{
				ColHi = Col;
				break;
				}
			}
		asserta(Sum >= RVA[ColX]);
		}

	unsigned ColXHi = ColX;
	for (unsigned Col = ColX + 1; Col < ColCount; ++Col)
		{
		char q = Q3Seq[Col];
		char a = A3Seq[Col];
		char b = B3Seq[Col];
		
		if (q == a && q == b && !isgap(q))
			ColXHi = Col;
		else
			break;
		}

	unsigned ColXLo = ColX;
	for (int Col = (int) ColX - 1; Col >= 0; --Col)
		{
		char q = Q3Seq[Col];
		char a = A3Seq[Col];
		char b = B3Seq[Col];
		
		if (q == a && q == b && !isgap(q))
			ColXLo = Col;
		else
			break;
		}

	unsigned IdQA = 0;
	unsigned IdQB = 0;
	unsigned IdAB = 0;
	unsigned NQA = 0;
	unsigned NQB = 0;
	unsigned NAB = 0;
	for (unsigned Col = 0; Col < ColCount; ++Col)
		{
		char q = Q3Seq[Col];
		char a = A3Seq[Col];
		char b = B3Seq[Col];

		if (!isgap(q) && !isgap(a))
			{
			++NQA;
			if (q == a)
				++IdQA;
			}

		if (!isgap(q) && !isgap(b))
			{
			++NQB;
			if (q == b)
				++IdQB;
			}

		if (!isgap(a) && !isgap(b))
			{
			++NAB;
			if (a == b)
				++IdAB;
			}
		}

	Hit.PctIdQA = Pct(IdQA, NQA);
	Hit.PctIdQB = Pct(IdQB, NQB);
	Hit.PctIdAB = Pct(IdAB, NAB);

	unsigned LIdQA = 0;
	unsigned LIdQB = 0;
	for (unsigned Col = ColLo; Col < ColXLo; ++Col)
		{
		char q = Q3Seq[Col];
		char a = A3Seq[Col];
		char b = B3Seq[Col];

		if (!isgap(q) && !isgap(a))
			{
			if (q == a)
				++LIdQA;
			}

		if (!isgap(q) && !isgap(b))
			{
			if (q == b)
				++LIdQB;
			}
		}

	unsigned RIdQA = 0;
	unsigned RIdQB = 0;
	for (unsigned Col = ColXHi+1; Col <= ColHi; ++Col)
		{
		char q = Q3Seq[Col];
		char a = A3Seq[Col];
		char b = B3Seq[Col];

		if (!isgap(q) && !isgap(a))
			{
			if (q == a)
				++RIdQA;
			}

		if (!isgap(q) && !isgap(b))
			{
			if (q == b)
				++RIdQB;
			}
		}

	unsigned IdDiffL = max(LIdQA, LIdQB) - min(LIdQA, LIdQB);
	unsigned IdDiffR = max(RIdQA, RIdQB) - min(RIdQA, RIdQB);
	unsigned MinIdDiff = min(IdDiffL, IdDiffR);
	unsigned ColRange = ColHi - ColLo + 1;
	if (opt_queryfract > 0.0f && float(ColRange)/float(QL) < opt_queryfract)
		return;

//	double Div = Pct(MinIdDiff, QSD.L);

#if	TRACE
	{
	Log("  Col  A Q B   ScoreA   ScoreB      LVA      LVB      RVA      RVB\n");
	Log("-----  - - -  -------  -------  -------  -------  -------  -------\n");
	for (unsigned Col = 0; Col < ColCount; ++Col)
		{
		if (ColScoresA[Col] == 0.0 && ColScoresB[Col] == 0.0)
			continue;

		char q = Q3Seq[Col];
		char a = A3Seq[Col];
		char b = B3Seq[Col];
		Log("%5u  %c %c %c", Col, a, q, b);

		if (ColScoresA[Col] == 0.0)
			Log("  %7.7s", "");
		else
			Log("  %7.1f", ColScoresA[Col]);

		if (ColScoresB[Col] == 0.0)
			Log("  %7.7s", "");
		else
			Log("  %7.1f", ColScoresB[Col]);

		Log("  %7.1f  %7.1f  %7.1f  %7.1f", LVA[Col], LVB[Col], RVA[Col], RVB[Col]);

		Log("\n");
		}
	Log("\n");
	Log("MaxSum %.1f, ColLo %u, ColXLo %u, ColX %u, ColXHi %u, ColHi %u, AF %c\n",
	  MaxSum, ColLo, ColXLo, ColX, ColXHi, ColHi, tof(FirstA));
	Log("  LIdQA %u, LIdQB %u, RIdQA %u, RIdQB %u\n", LIdQA, LIdQB, RIdQA, RIdQB);
	}
#endif

	string Q3L;
	string A3L;
	string B3L;
	for (unsigned Col = ColLo; Col <= ColHi; ++Col)
		{
		char q = Q3[Col];
		char a = A3[Col];
		char b = B3[Col];

		Q3L += q;
		A3L += a;
		B3L += b;
		}

	AlignChimeGlobal3(Q3L, A3L, B3L, QLabel, ALabel, BLabel, Hit);

#if	0
// CS SNPs
	Hit.CS_LY = 0;
	Hit.CS_LN = 0;
	Hit.CS_RY = 0;
	Hit.CS_RN = 0;
	Hit.CS_LA = 0;
	Hit.CS_RA = 0;
	for (unsigned Col = ColLo; Col <= ColHi; ++Col)
		{
		char q = Q3Seq[Col];
		char a = A3Seq[Col];
		char b = B3Seq[Col];
		if (q == a && q == b && a == b)
			continue;
		if (isgap(q) || isgap(a) || isgap(b))
			continue;
		if (Col > 0 && (isgap(Q3Seq[Col-1]) || isgap(A3Seq[Col-1]) || isgap(B3Seq[Col-1])))
			continue;
		if (Col + 1 < ColCount && (isgap(Q3Seq[Col+1]) || isgap(A3Seq[Col+1]) || isgap(B3Seq[Col+1])))
			continue;

		if (!FirstA)
			swap(a, b);

		if (Col < ColXLo)
			{
			if (q == a && q != b)
				++Hit.CS_LY;
			else if (q == b && q != a)
				++Hit.CS_LN;
			else
				++Hit.CS_LA;
			}
		else if (Col > ColXHi)
			{
			if (q == b && q != a)
				++Hit.CS_RY;
			else if (q == a && q != b)
				++Hit.CS_RN;
			else
				++Hit.CS_RA;
			}
		}

	double ScoreL = GetScore2(Hit.CS_LY, Hit.CS_LN, Hit.CS_LA);
	double ScoreR = GetScore2(Hit.CS_RY, Hit.CS_RN, Hit.CS_RA);
	Hit.Score = ScoreL*ScoreR;

	//Hit.QSD = QSD;
	//if (FirstA)
	//	{
	//	Hit.ASD = ASD;
	//	Hit.BSD = BSD;
	//	Hit.PathQA = PathQA;
	//	Hit.PathQB = PathQB;
	//	}
	//else
	//	{
	//	Hit.ASD = BSD;
	//	Hit.BSD = ASD;
	//	}

	//Hit.ColLo = ColLo;
	//Hit.ColXLo = ColXLo;
	//Hit.ColXHi = ColXHi;
	//Hit.ColHi = ColHi;
	//Hit.Div = Div;

//	Hit.LogMe();
#endif
	}
