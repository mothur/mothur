#include "myutils.h"
#include <float.h>	// for FLT_MAX
#include "mx.h"
#include "alnparams.h"
#include "hsp.h"

#define TEST	0

void SetBLOSUM62();
void SetNucSubstMx(double Match, double Mismatch);
void ReadSubstMx(const string &FileName, Mx<float> &Mxf);

extern Mx<float> g_SubstMxf;
extern float **g_SubstMx;

void AlnParams::Clear()
	{
	SubstMxName = 0;
	LocalOpen = OBVIOUSLY_WRONG_PENALTY;
	LocalExt = OBVIOUSLY_WRONG_PENALTY;
	OpenA = OBVIOUSLY_WRONG_PENALTY;
	OpenB = OBVIOUSLY_WRONG_PENALTY;
	ExtA = OBVIOUSLY_WRONG_PENALTY;
	ExtB = OBVIOUSLY_WRONG_PENALTY;
	LOpenA = OBVIOUSLY_WRONG_PENALTY;
	LOpenB = OBVIOUSLY_WRONG_PENALTY;
	ROpenA = OBVIOUSLY_WRONG_PENALTY;
	ROpenB = OBVIOUSLY_WRONG_PENALTY;
	LExtA = OBVIOUSLY_WRONG_PENALTY;
	LExtB = OBVIOUSLY_WRONG_PENALTY;
	RExtA = OBVIOUSLY_WRONG_PENALTY;
	RExtB = OBVIOUSLY_WRONG_PENALTY;
	Nucleo = false;
	NucleoSet = false;
	}

bool AlnParams::Is2() const
	{
	float g = OpenA;
	float e = ExtA;
	if (OpenB != g || LOpenA != g || LOpenB != g || ROpenA != g || ROpenB != g)
		return false;
	if (ExtB != e || LExtA != e || LExtB != e || RExtA != e || RExtB != e)
		return false;
	return true;
	}

bool AlnParams::Is4() const
	{
	float g = OpenA;
	float tg = LOpenA;
	float e = ExtA;
	float te = LExtA;
	if (OpenB != g || LOpenA != tg || LOpenB != tg || ROpenA != tg || ROpenB != tg)
		return false;
	if (ExtB != e || LExtA != te || LExtB != te || RExtA != te || RExtB != te)
		return false;
	return true;
	}

const char *AlnParams::GetType() const
	{
	if (Is2())
		return "2";
	else if (Is4())
		return "4";
	return "12";
	}

void AlnParams::Init2(const float * const *Mx, float Open, float Ext)
	{
	SubstMx = Mx;
	OpenA = OpenB = LOpenA = LOpenB = ROpenA = ROpenB = Open;
	ExtA = ExtB = LExtA = LExtB = RExtA = RExtB = Ext;
	}

void AlnParams::SetLocal(float Open, float Ext)
	{
	LocalOpen = Open;
	LocalExt = Ext;
	}

void AlnParams::Init4(const float * const *Mx, float Open, float Ext,
  float TermOpen, float TermExt)
	{
	SubstMx = Mx;
	OpenA = OpenB = Open;
	LOpenA = LOpenB = ROpenA = ROpenB = TermOpen;
	ExtA = ExtB = Ext;
	LExtA = LExtB = RExtA = RExtB = TermExt;
	}

void AlnParams::Init(const AlnParams &AP, const HSPData &HSP,
  unsigned LA, unsigned LB)
	{
	SubstMx = AP.SubstMx;
	OpenA = AP.OpenA;
	OpenB = AP.OpenB;
	ExtA = AP.ExtA;
	ExtB = AP.ExtB;

	if (HSP.LeftA())
		{
		LOpenA = AP.LOpenA;
		LExtA = AP.LExtA;
		}
	else
		{
		LOpenA = AP.OpenA;
		LExtA = AP.ExtA;
		}

	if (HSP.LeftB())
		{
		LOpenB = AP.LOpenB;
		LExtB = AP.LExtB;
		}
	else
		{
		LOpenB = AP.OpenB;
		LExtB = AP.ExtB;
		}

	if (HSP.RightA(LA))
		{
		ROpenA = AP.ROpenA;
		RExtA = AP.RExtA;
		}
	else
		{
		ROpenA = AP.OpenA;
		RExtA = AP.ExtA;
		}

	if (HSP.RightB(LB))
		{
		ROpenB = AP.ROpenB;
		RExtB = AP.RExtB;
		}
	else
		{
		ROpenB = AP.OpenB;
		RExtB = AP.ExtB;
		}
	}

void AlnParams::LogMe() const
	{
	Log("AlnParams(%s)", GetType());
	if (Is2())
		Log(" g=%.1f e=%.1f", -OpenA, -ExtA);
	else if (Is4())
		Log(" g=%.1f tg=%.1f e=%.1f te=%.1f", -OpenA, -ExtA, -LOpenA, -LExtA);
	else
		Log(
" gA=%.1f gB=%.1f gAL=%.1f gBL=%.1f gAR=%.1f gBR=%.1f eA=%.1f eB=%.1f eAL=%.1f eBL=%.1f eAR=%.1f eBR=%.1f",
		  OpenA, OpenB, LOpenA, LOpenB, ROpenA, ROpenB, ExtA, ExtB, LExtA, LExtB, RExtA, RExtB);
	Log("\n");
	}

/***
Open/Ext format string is one or more:
	[<flag><flag>...]<value>

Value is (positive) penalty or * (disabled).
Flag is:
	Q		Query.
	T		Target sequence.
	I		Internal gaps (defafault internal and terminal).
	E		End gaps (default internal and terminal).
	L		Left end.
	R		Right end.
***/

static void ParseGapStr(const string &s,
  float &QI, float &QL, float &QR,
  float &TI, float &TL, float &TR)
	{
	if (s.empty())
		return;

	bool Q = false;
	bool T = false;
	bool I = false;
	bool E = false;
	bool L = false;
	bool R = false;

	const unsigned K = SIZE(s);
	unsigned Dec = 0;
	float Value = FLT_MAX;
	for (unsigned i = 0; i <= K; ++i)
		{
		char c = s.c_str()[i];
		if (c == 0 || c == '/')
			{
			if (Value == FLT_MAX)
				Die("Invalid gap penalty string, missing penalty '%s'", s.c_str());
			if (!Q && !T && !I && !E && !L && !R)
				{
				Q = true;
				T = true;
				L = true;
				R = true;
				I = true;
				}

			if (!E && !I && !L && !R)
				{
				E = false;
				I = true;
				L = true;
				R = true;
				}

			if (E)
				{
				if (L || R)
					Die("Invalid gap penalty string (E and L or R) '%s'", s.c_str());
				L = true;
				R = true;
				}

			if (!Q && !T)
				{
				Q = true;
				T = true;
				}

			if (Q && L)
				QL = -Value;
			if (Q && R)
				QR = -Value;
			if (Q && I)
				QI = -Value;
			if (T && L)
				TL = -Value;
			if (T && R)
				TR = -Value;
			if (T && I)
				TI = -Value;
			
			Value = FLT_MAX;
			Dec = 0;
			Q = false;
			T = false;
			I = false;
			E = false;
			L = false;
			R = false;
			}
		else if (c == '*')
			{
			if (Value != FLT_MAX)
				Die("Invalid gap penalty (* in floating point number) '%s'", s.c_str());
			Value = -MINUS_INFINITY;
			}
		else if (isdigit(c))
			{
			if (Value == -MINUS_INFINITY)
				Die("Invalid gap penalty (* in floating point number) '%s'", s.c_str());
			if (Value == FLT_MAX)
				Value = 0.0;
			if (Dec > 0)
				{
				Dec *= 10;
				Value += float(c - '0')/Dec;
				}
			else
				Value = Value*10 + (c - '0');
			}
		else if (c == '.')
			{
			if (Dec > 0)
				Die("Invalid gap penalty (two decimal points) '%s'", s.c_str());
			Dec = 1;
			}
		else
			{
			switch (c)
				{
			case 'Q':
				Q = true;
				break;
			case 'T':
				T = true;
				break;
			case 'I':
				I = true;
				break;
			case 'L':
				L = true;
				break;
			case 'R':
				R = true;
				break;
			case 'E':
				E = true;
				break;
			default:
				Die("Invalid char '%c' in gap penalty string '%s'", c, s.c_str());
				}
			}
		}
	}

void AlnParams::SetPenalties(const string &OpenStr, const string &ExtStr)
	{
	ParseGapStr(OpenStr, OpenA, LOpenA, ROpenA, OpenB, LOpenB, ROpenB);
	ParseGapStr(ExtStr, ExtA, LExtA, RExtA, ExtB, LExtB, RExtB);
	}

void AlnParams::SetMxFromCmdLine(bool IsNucleo)
	{
	if (IsNucleo)
		SetNucSubstMx(opt_match, opt_mismatch);
	else
		{
		if (opt_matrix == "")
			{
			SubstMxName = "BLOSUM62";
			SetBLOSUM62();
			}
		else
			{
			ReadSubstMx(opt_matrix, g_SubstMxf);
			g_SubstMx = g_SubstMxf.GetData();
			g_SubstMxf.LogMe();
			SubstMxName = opt_matrix.c_str();
			}
		}
	SubstMx = g_SubstMx;
	asserta(SubstMx != 0);
	}

void AlnParams::InitFromCmdLine(bool IsNucleo)
	{
	Clear();
	Nucleo = IsNucleo;
	NucleoSet = true;

	SetMxFromCmdLine(IsNucleo);

// Local
	if (optset_lopen || optset_lext)
		{
		if (!optset_lopen || !optset_lext)
			Die("Must set both --lopen and --lext");
		if (opt_lopen < 0.0 || opt_lext < 0.0)
			Die("Invalid --lopen/--lext, gap penalties must be >= 0");
		SetLocal(float(-opt_lopen), float(-opt_lext));
		}
	else
		{
	// Same penalties, if-statement to note could differ.
		if (IsNucleo)
			SetLocal(-10.0f, -1.0f);
		else
			SetLocal(-10.0f, -1.0f);
		}

// Global
	if (IsNucleo)
		Init4(g_SubstMx, -10.0, -1.0, -0.5, -0.5);
	else
		Init4(g_SubstMx, -17.0, -1.0, -0.5, -0.5);
	SetPenalties(opt_gapopen, opt_gapext);
	}

float AlnParams::GetLocalOpen() const
	{
	return LocalOpen;
	}

float AlnParams::GetLocalExt() const
	{
	return LocalExt;
	}

bool AlnParams::GetIsNucleo() const
	{
	asserta(NucleoSet);
	return Nucleo;
	}

unsigned GetWindexWordLength(bool Nucleo)
	{
	if (optset_w)
		return opt_w;

	if (Nucleo)
		return 8;
	else
		return 5;
	}

#if	TEST
static void Test1(const string &os, const string &es)
	{
	AlnParams AP;
	Log("\n");
	Log("OpenStr %s\n", os.c_str());
	Log(" ExtStr %s\n", es.c_str());
	AP.SetPenalties(os, es);
	AP.LogMe();
	}

void TestGapStr()
	{
	Test1("17I/0.5E", "1I/0.5E");
	Test1("17I/0.5L/0.4R", "1Q/2T");
	Test1("1QL/2QR/3QI/4TL/5TR/6TI", ".1QL/.2QR/.3QI/.4TL/.5TR/.6TI");
	}
#endif // TEST
