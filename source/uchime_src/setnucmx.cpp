#include "myutils.h"
#include "mx.h"

Mx<float> g_SubstMxf;
float **g_SubstMx;

static const char Alphabet[] = "ACGTU";

void SetNucSubstMx(double Match, double Mismatch)
	{
	static bool Done = false;
	if (Done)
		return;
	Done = true;

	if (Match <= 0.0)
		Die("Match score should be +ve");
	if (Mismatch >= 0.0)
		Die("Mismatch score should be -ve");

	unsigned N = unsigned(strlen(Alphabet));

	g_SubstMxf.Alloc("NUCMX", 256, 256);
	strcpy(g_SubstMxf.m_Alpha, "ACGT");
	g_SubstMxf.Init(0);
	g_SubstMx = g_SubstMxf.GetData();
	for (unsigned i = 0; i < N; ++i)
		{
		for (unsigned j = 0; j < N; ++j)
			{
			float v = float(i == j ? Match : Mismatch);

			byte ui = (byte) toupper(Alphabet[i]);
			byte uj = (byte) toupper(Alphabet[j]);
			byte li = (byte) tolower(ui);
			byte lj = (byte) tolower(uj);
			ui = (byte) toupper(ui);
			uj = (byte) toupper(uj);

			g_SubstMx[ui][uj] = v;
			g_SubstMx[uj][ui] = v;

			g_SubstMx[ui][lj] = v;
			g_SubstMx[uj][li] = v;

			g_SubstMx[li][uj] = v;
			g_SubstMx[lj][ui] = v;

			g_SubstMx[li][lj] = v;
			g_SubstMx[lj][li] = v;
			}
		}

	for (unsigned j = 0; j < N; ++j)
		{
		float v = 0.0f;

		byte ui = (byte) 'N';
		byte uj = (byte) toupper(Alphabet[j]);
		byte li = (byte) 'n';
		byte lj = (byte) tolower(uj);
		ui = (byte) toupper(ui);
		uj = (byte) toupper(uj);

		g_SubstMx[ui][uj] = v;
		g_SubstMx[uj][ui] = v;

		g_SubstMx[ui][lj] = v;
		g_SubstMx[uj][li] = v;

		g_SubstMx[li][uj] = v;
		g_SubstMx[lj][ui] = v;

		g_SubstMx[li][lj] = v;
		g_SubstMx[lj][li] = v;
		}
	}
