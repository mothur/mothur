//uchime by Robert C. Edgar http://drive5.com/uchime This code is donated to the public domain.

#include "myutils.h"
#include "seq.h"
#include "chime.h"
#include "dp.h"

#define TRACE		0
#define TRACE_BS	0

void Make3Way(const SeqData &SDQ, const SeqData &SDA, const SeqData &SDB,
  const string &PathQA, const string &PathQB,
  string &Q3, string &A3, string &B3);

void AlignChimeLocal3(const string &Q3, const string &A3, const string &B3,
  const string &QLabel, const string &ALabel, const string &BLabel,
  ChimeHit2 &Hit);

double GetScore2(double Y, double N, double A)
	{
	return Y/(opt_xn*(N + opt_dn) + opt_xa*A);
	}

void AlignChimeGlobal3(const string &Q3, const string &A3, const string &B3,
  const string &QLabel, const string &ALabel, const string &BLabel,
  ChimeHit2 &Hit)
	{
	Hit.Clear();
	Hit.QLabel = QLabel;

	const byte *Q3Seq = (const byte *) Q3.c_str();
	const byte *A3Seq = (const byte *) A3.c_str();
	const byte *B3Seq = (const byte *) B3.c_str();

	const unsigned ColCount = SIZE(Q3);
	asserta(SIZE(A3) == ColCount && SIZE(B3) == ColCount);

#if	TRACE
	Log("Q %5u %*.*s\n", ColCount, ColCount, ColCount, Q3Seq);
	Log("A %5u %*.*s\n", ColCount, ColCount, ColCount, A3Seq);
	Log("B %5u %*.*s\n", ColCount, ColCount, ColCount, B3Seq);
#endif

// Discard terminal gaps
	unsigned ColLo = UINT_MAX;
	unsigned ColHi = UINT_MAX;
	for (unsigned Col = 2; Col + 2 < ColCount; ++Col)
		{
		char q = Q3Seq[Col];
		char a = A3Seq[Col];
		char b = B3Seq[Col];

		if (isacgt(q) && isacgt(a) && isacgt(b))
			{
			if (ColLo == UINT_MAX)
				ColLo = Col;
			ColHi = Col;
			}
		}

	if (ColLo == UINT_MAX)
		return;

	unsigned QPos = 0;
	unsigned APos = 0;
	unsigned BPos = 0;
	unsigned DiffCount = 0;

	vector<unsigned> ColToQPos(ColLo, UINT_MAX);
	vector<unsigned> AccumCount(ColLo, UINT_MAX);
	vector<unsigned> AccumSameA(ColLo, UINT_MAX);
	vector<unsigned> AccumSameB(ColLo, UINT_MAX);
	vector<unsigned> AccumForA(ColLo, UINT_MAX);
	vector<unsigned> AccumForB(ColLo, UINT_MAX);
	vector<unsigned> AccumAbstain(ColLo, UINT_MAX);
	vector<unsigned> AccumAgainst(ColLo, UINT_MAX);

	unsigned SumSameA = 0;
	unsigned SumSameB = 0;
	unsigned SumSameAB = 0;
	unsigned Sum = 0;
	unsigned SumForA = 0;
	unsigned SumForB = 0;
	unsigned SumAbstain = 0;
	unsigned SumAgainst = 0;
	for (unsigned Col = ColLo; Col <= ColHi; ++Col)
		{
		char q = Q3Seq[Col];
		char a = A3Seq[Col];
		char b = B3Seq[Col];

		if (isacgt(q) && isacgt(a) && isacgt(b))
			{
			if (q == a)
				++SumSameA;
			if (q == b)
				++SumSameB;
			if (a == b)
				++SumSameAB;
			if (q == a && q != b)
				++SumForA;
			if (q == b && q != a)
				++SumForB;
			if (a == b && q != a)
				++SumAgainst;
			if (q != a && q != b)
				++SumAbstain;
			++Sum;
			}

		ColToQPos.push_back(QPos);
		AccumSameA.push_back(SumSameA);
		AccumSameB.push_back(SumSameB);
		AccumCount.push_back(Sum);
		AccumForA.push_back(SumForA);
		AccumForB.push_back(SumForB);
		AccumAbstain.push_back(SumAbstain);
		AccumAgainst.push_back(SumAgainst);

		if (q != '-')
			++QPos;
		if (a != '-')
			++APos;
		if (b != '-')
			++BPos;
		}

	asserta(SIZE(ColToQPos) == ColHi+1);
	asserta(SIZE(AccumSameA) == ColHi+1);
	asserta(SIZE(AccumSameB) == ColHi+1);
	asserta(SIZE(AccumAbstain) == ColHi+1);
	asserta(SIZE(AccumAgainst) == ColHi+1);

	double IdQA = double(SumSameA)/Sum;
	double IdQB = double(SumSameB)/Sum;
	double IdAB = double(SumSameAB)/Sum;
	double MaxId = max(IdQA, IdQB);

#if	TRACE
	Log("IdQA=%.1f%% IdQB=%.1f%% IdAB=%.1f\n", IdQA*100.0, IdQB*100.0, IdAB*100.0);
	Log("\n");
	Log("    x  AQB   IdAL   IdBL   IdAR   IdBR   DivAB   DivBA    YAL    YBL    YAR    YBR    AbL    AbR  ScoreAB  ScoreAB    XLo    Xhi\n");
	Log("-----  ---  -----  -----  -----  -----  ------  ------  -----  -----  -----  -----  -----  -----  -------  -------  -----  -----\n");
#endif
	unsigned BestXLo = UINT_MAX;
	unsigned BestXHi = UINT_MAX;
	double BestDiv = 0.0;
	double BestIdQM = 0.0;
	double BestScore = 0.0;

// Find range of cols BestXLo..BestXHi that maximizes score
	bool FirstA = false;

// NOTE: Must be < ColHi not <= because use Col+1 below
	for (unsigned Col = ColLo; Col < ColHi; ++Col)
		{
		char q = Q3Seq[Col];
		char a = A3Seq[Col];
		char b = B3Seq[Col];

		unsigned SameAL = AccumSameA[Col];
		unsigned SameBL = AccumSameB[Col];
		unsigned SameAR = SumSameA - AccumSameA[Col];
		unsigned SameBR = SumSameB - AccumSameB[Col];

		double IdAB = double(SameAL + SameBR)/Sum;
		double IdBA = double(SameBL + SameAR)/Sum;

		unsigned ForAL = AccumForA[Col];
		unsigned ForBL = AccumForB[Col];
		unsigned ForAR = SumForA - AccumForA[Col+1];
		unsigned ForBR = SumForB - AccumForB[Col+1];
		unsigned AbL = AccumAbstain[Col];
		unsigned AbR = SumAbstain - AccumAbstain[Col+1];

		double ScoreAB = GetScore2(ForAL, ForBL, AbL)*GetScore2(ForBR, ForAR, AbR);
		double ScoreBA = GetScore2(ForBL, ForAL, AbL)*GetScore2(ForAR, ForBR, AbR);
	
		double DivAB = IdAB/MaxId;
		double DivBA = IdBA/MaxId;
		double MaxDiv = max(DivAB, DivBA);

		//if (MaxDiv > BestDiv)
		//	{
		//	BestDiv = MaxDiv;
		//	BestXLo = Col;
		//	BestXHi = Col;
		//	FirstA = (DivAB > DivBA);
		//	if (FirstA)
		//		BestIdQM = IdAB;
		//	else
		//		BestIdQM = IdBA;
		//	}
		//else if (MaxDiv == BestDiv)
		//	BestXHi = Col;

		double MaxScore = max(ScoreAB, ScoreBA);
		if (MaxScore > BestScore)
			{
			BestScore = MaxScore;
			BestXLo = Col;
			BestXHi = Col;
			FirstA = (ScoreAB > ScoreBA);
			if (FirstA)
				BestIdQM = IdAB;
			else
				BestIdQM = IdBA;
			if (MaxDiv > BestDiv)
				BestDiv = MaxDiv;
			}
		else if (MaxScore == BestScore)
			{
			BestXHi = Col;
			if (MaxDiv > BestDiv)
				BestDiv = MaxDiv;
			}

#if	TRACE
		{
		Log("%5u", Col);
		char q = Q3Seq[Col];
		char a = A3Seq[Col];
		char b = B3Seq[Col];
		Log("  %c%c%c", a, q, b);
		Log("  %5u", SameAL);
		Log("  %5u", SameBL);
		Log("  %5u", SameAR);
		Log("  %5u", SameBR);
		Log("  %5.4f", DivAB);
		Log("  %5.4f", DivBA);
		Log("  %5u", ForAL);
		Log("  %5u", ForBL);
		Log("  %5u", ForAR);
		Log("  %5u", ForBR);
		Log("  %5u", AbL);
		Log("  %5u", AbR);
		Log("  %7.4f", ScoreAB);
		Log("  %7.4f", ScoreBA);
		if (BestXLo != UINT_MAX)
			Log("  %5u", BestXLo);
		if (BestXHi != UINT_MAX)
			Log("  %5u", BestXHi);
		Log("\n");
		}
#endif
		}

	if (BestXLo == UINT_MAX)
		{
#if	TRACE
		Log("\n");
		Log("No crossover found.\n");
#endif
		return;
		}
#if	TRACE
	Log("BestX col %u - %u\n", BestXLo, BestXHi);
#endif

// Find maximum region of identity within BestXLo..BestXHi
	unsigned ColXLo = (BestXLo + BestXHi)/2;
	unsigned ColXHi = ColXLo;
	unsigned SegLo = UINT_MAX;
	unsigned SegHi = UINT_MAX;
	for (unsigned Col = BestXLo; Col <= BestXHi; ++Col)
		{
		char q = Q3Seq[Col];
		char a = A3Seq[Col];
		char b = B3Seq[Col];

		if (q == a && q == b)
			{
			if (SegLo == UINT_MAX)
				SegLo = Col;
			SegHi = Col;
			}
		else
			{
			unsigned SegLength = SegHi - SegLo + 1;
			unsigned BestSegLength = ColXHi - ColXLo + 1;
			if (SegLength > BestSegLength)
				{
				ColXLo = SegLo;
				ColXHi = SegHi;
				}
			SegLo = UINT_MAX;
			SegHi = UINT_MAX;
			}
		}
	unsigned SegLength = SegHi - SegLo + 1;
	unsigned BestSegLength = ColXHi - ColXLo + 1;
	if (SegLength > BestSegLength)
		{
		ColXLo = SegLo;
		ColXHi = SegHi;
		}

	QPos = 0;
	for (unsigned x = 0; x < ColCount; ++x)
		{
		if (x == ColXLo)
			Hit.QXLo = QPos;
		else if (x == ColXHi)
			{
			Hit.QXHi = QPos;
			break;
			}
		char q = Q3Seq[x];
		if (q != '-')
			++QPos;
		}

	Hit.ColXLo = ColXLo;
	Hit.ColXHi = ColXHi;

	//if (FirstA)
	//	{
	//	Hit.LY = AccumForA[ColXLo];
	//	Hit.LN = AccumForB[ColXLo];

	//	Hit.RY = SumForB - AccumForB[ColXHi];
	//	Hit.RN = SumForA - AccumForA[ColXHi];
	//	}
	//else
	//	{
	//	Hit.LY = AccumForB[ColXLo];
	//	Hit.LN = AccumForA[ColXLo];
	//	Hit.RY = SumForA - AccumForA[ColXHi];
	//	Hit.RN = SumForB - AccumForB[ColXHi];
	//	}

	//Hit.LA = AccumAgainst[ColXLo];
	//Hit.LD = AccumAbstain[ColXLo];

	//Hit.RA = SumAgainst - AccumAgainst[ColXHi];
	//Hit.RD = SumAbstain - AccumAbstain[ColXHi];

	Hit.PctIdAB = IdAB*100.0;
	Hit.PctIdQM = BestIdQM*100.0;

	Hit.Div = (BestDiv - 1.0)*100.0;

	//Hit.QSD = QSD;
	Hit.Q3 = Q3;
	Hit.QLabel = QLabel;
	if (FirstA)
		{
		//Hit.ASD = ASD;
		//Hit.BSD = BSD;
		//Hit.PathQA = PathQA;
		//Hit.PathQB = PathQB;
		Hit.A3 = A3;
		Hit.B3 = B3;
		Hit.ALabel = ALabel;
		Hit.BLabel = BLabel;
		Hit.PctIdQA = IdQA*100.0;
		Hit.PctIdQB = IdQB*100.0;
		}
	else
		{
		Hit.A3 = B3;
		Hit.B3 = A3;
		Hit.ALabel = BLabel;
		Hit.BLabel = ALabel;
		Hit.PctIdQA = IdQB*100.0;
		Hit.PctIdQB = IdQA*100.0;
		}

// CS SNPs
	Hit.CS_LY = 0;
	Hit.CS_LN = 0;
	Hit.CS_RY = 0;
	Hit.CS_RN = 0;
	Hit.CS_LA = 0;
	Hit.CS_RA = 0;

	//vector<float> Cons;
	//for (unsigned Col = 0; Col < ColCount; ++Col)
	//	{
	//	char q = Q3Seq[Col];
	//	char a = A3Seq[Col];
	//	char b = B3Seq[Col];
	//	if (q == a && q == b && a == b)
	//		{
	//		Cons.push_back(1.0f);
	//		continue;
	//		}

	//	bool gapq = isgap(q);
	//	bool gapa = isgap(a);
	//	bool gapb = isgap(b);

	//	if (!gapq && !gapa && !gapb)
	//		{
	//		if (q == a || q == b || a == b)
	//			Cons.push_back(0.75);
	//		else
	//			Cons.push_back(0.5);
	//		}
	//	else
	//		{
	//		if (!gapa && (a == b || a == q))
	//			Cons.push_back(0.5f);
	//		else if (!gapb && b == q)
	//			Cons.push_back(0.5f);
	//		else
	//			Cons.push_back(0.0f);
	//		}
	//	}

	//float fLY = 0.0f;
	//float fLN = 0.0f;
	//float fLA = 0.0f;
	//float fRY = 0.0f;
	//float fRN = 0.0f;
	//float fRA = 0.0f;
	for (unsigned Col = ColLo; Col <= ColHi; ++Col)
		{
		char q = Q3Seq[Col];
		char a = A3Seq[Col];
		char b = B3Seq[Col];
		if (q == a && q == b && a == b)
			continue;

		unsigned ngaps = 0;
		if (isgap(q))
			++ngaps;
		if (isgap(a))
			++ngaps;
		if (isgap(b))
			++ngaps;

		if (opt_skipgaps)
			{
			if (ngaps == 3)
				continue;
			}
		else
			{
			if (ngaps == 2)
				continue;
			}

		if (!FirstA)
			swap(a, b);

		//float AvgCons = (Cons[Col-2] + Cons[Col-1] + Cons[Col+1] + Cons[Col+2])/4;
		//if (Col < ColXLo)
		//	{
		//	if (q == a && q != b)
		//		fLY += AvgCons;
		//	else if (q == b && q != a)
		//		fLN += AvgCons;
		//	else
		//		fLA += AvgCons;
		//	}
		//else if (Col > ColXHi)
		//	{
		//	if (q == b && q != a)
		//		fRY += AvgCons;
		//	else if (q == a && q != b)
		//		fRN += AvgCons;
		//	else
		//		fRA += AvgCons;
		//	}

		if (opt_skipgaps2)
			{
			if (Col > 0 && (isgap(Q3Seq[Col-1]) || isgap(A3Seq[Col-1]) || isgap(B3Seq[Col-1])))
				continue;
			if (Col + 1 < ColCount && (isgap(Q3Seq[Col+1]) || isgap(A3Seq[Col+1]) || isgap(B3Seq[Col+1])))
				continue;
			}

		//if (Col > 0 && isgap(Q3Seq[Col-1]))
			//continue;
		//if (Col + 1 < ColCount && isgap(Q3Seq[Col+1]))
		//	continue;

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

	extern bool g_UchimeDeNovo;

	//if (0)//g_UchimeDeNovo)
	//	{
	//	double AbQ = GetAbFromLabel(QLabel.c_str());
	//	double AbA = GetAbFromLabel(ALabel.c_str());
	//	double AbB = GetAbFromLabel(BLabel.c_str());
	//	if (AbQ > 0.0 && AbA > 0.0 && AbB > 0.0)
	//		{
	//		double MinAb = min(AbA, AbB);
	//		double Ratio = MinAb/AbQ;
	//		double t = Ratio - opt_abx;
	//	//	double Factor = 2.0/(1.0 + exp(-t));
	//		double Factor = min(Ratio, opt_abx)/opt_abx;
	//		if (opt_verbose)
	//			Log("Score %.4f Ab factor %.4f >%s\n", Hit.Score, Factor, QLabel.c_str());
	//		Hit.Score *= Factor;
	//		}
	//	}

	extern FILE *g_fUChimeAlns;
	if (g_fUChimeAlns != 0 && Hit.Div > 0.0)
		{
		void WriteChimeHitX(FILE *f, const ChimeHit2 &Hit);
		WriteChimeHitX(g_fUChimeAlns, Hit);
		}
	}

void AlignChime3(const string &Q3, const string &A3, const string &B3,
  const string &QLabel, const string &ALabel, const string &BLabel,
  ChimeHit2 &Hit)
	{
	if (opt_ucl)
		AlignChimeLocal3(Q3, A3, B3, QLabel, ALabel, BLabel, Hit);
	else
		AlignChimeGlobal3(Q3, A3, B3, QLabel, ALabel, BLabel, Hit);
	}

static void StripGaps(const byte *Seq, unsigned L, string &s)
	{
	s.clear();
	for (unsigned i = 0; i < L; ++i)
		{
		char c = Seq[i];
		if (!isgap(c))
			s.push_back(c);
		}
	}

static void StripGapsAlloc(const SeqData &SDIn, SeqData &SDOut)
	{
	SDOut = SDIn;
	byte *s = myalloc(byte, SDIn.L);
	unsigned k = 0;
	for (unsigned i = 0; i < SDIn.L; ++i)
		{
		char c = SDIn.Seq[i];
		if (!isgap(c))
			s[k++] = toupper(c);
		}
	SDOut.Seq = s;
	SDOut.L = k;
	}

void AlignChime(const SeqData &QSD, const SeqData &ASD, const SeqData &BSD,
  const string &PathQA, const string &PathQB, ChimeHit2 &Hit)
	{
	//if (opt_ucl)
	//	{
	//	AlignChimeLocal(QSD, ASD, BSD, PathQA, PathQB, Hit);
	//	return;
	//	}

	string Q3;
	string A3;
	string B3;
	Make3Way(QSD, ASD, BSD, PathQA, PathQB, Q3, A3, B3);

	AlignChime3(Q3, A3, B3, QSD.Label, ASD.Label, BSD.Label, Hit);
	}

void AlignChime3SDRealign(const SeqData &QSD3, const SeqData &ASD3, const SeqData &BSD3,
  ChimeHit2 &Hit)
	{
	SeqData QSD;
	SeqData ASD;
	SeqData BSD;
	StripGapsAlloc(QSD3, QSD);
	StripGapsAlloc(ASD3, ASD);
	StripGapsAlloc(BSD3, BSD);

	string PathQA;
	string PathQB;
	bool FoundQA = GlobalAlign(QSD, ASD, PathQA);
	bool FoundQB = GlobalAlign(QSD, BSD, PathQB);
	if (!FoundQA || !FoundQB)
		{
		Hit.Clear();
		Hit.QLabel = QSD3.Label;
		return;
		}

	AlignChime(QSD, ASD, BSD, PathQA, PathQB, Hit);

	myfree((void *) QSD.Seq);
	myfree((void *) ASD.Seq);
	myfree((void *) BSD.Seq);
	}

void AlignChime3SD(const SeqData &QSD3, const SeqData &ASD3, const SeqData &BSD3,
  ChimeHit2 &Hit)
	{
	if (opt_realign)
		{
		AlignChime3SDRealign(QSD3, ASD3, BSD3, Hit);
		return;
		}

	string Q3;
	string A3;
	string B3;

	const unsigned ColCount = QSD3.L;
	asserta(ASD3.L == ColCount && BSD3.L == ColCount);

	Q3.reserve(ColCount);
	A3.reserve(ColCount);
	B3.reserve(ColCount);

	const byte *QS = QSD3.Seq;
	const byte *AS = ASD3.Seq;
	const byte *BS = BSD3.Seq;
	for (unsigned Col = 0; Col < ColCount; ++Col)
		{
		byte q = toupper(QS[Col]);
		byte a = toupper(AS[Col]);
		byte b = toupper(BS[Col]);

		if (isgap(q) && isgap(a) && isgap(b))
			continue;

		Q3.push_back(q);
		A3.push_back(a);
		B3.push_back(b);
		}

	AlignChime3(Q3, A3, B3, QSD3.Label, ASD3.Label, BSD3.Label, Hit);
	}
