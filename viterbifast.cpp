//uchime by Robert C. Edgar http://drive5.com/uchime This code is donated to the public domain.

#include "dp.h"
#include "out.h"
#include "evalue.h"

#define CMP_SIMPLE	0

#if	SAVE_FAST
static Mx<float> g_MxDPM;
static Mx<float> g_MxDPD;
static Mx<float> g_MxDPI;

static Mx<char> g_MxTBM;
static Mx<char> g_MxTBD;
static Mx<char> g_MxTBI;

static float **g_DPM;
static float **g_DPD;
static float **g_DPI;

static char **g_TBM;
static char **g_TBD;
static char **g_TBI;

#if	CMP_SIMPLE
static Mx<float> *g_DPMSimpleMx;
static Mx<float> *g_DPDSimpleMx;
static Mx<float> *g_DPISimpleMx;
static float **g_DPMSimple;
static float **g_DPDSimple;
static float **g_DPISimple;

#define cmpm(i, j, x)	{ if (!feq(x, g_DPMSimple[i][j])) \
							{ \
							Die("%s:%d %.1f != DPMSimple[%u][%u] = %.1f", \
							  __FILE__, __LINE__, x, i, j, g_DPMSimple[i][j]); \
							} \
						}

#define cmpd(i, j, x)	{ if (!feq(x, g_DPDSimple[i][j])) \
							{ \
							Die("%s:%d %.1f != DPMSimple[%u][%u] = %.1f", \
							  __FILE__, __LINE__, x, i, j, g_DPDSimple[i][j]); \
							} \
						}

#define cmpi(i, j, x)	{ if (!feq(x, g_DPISimple[i][j])) \
							{ \
							Die("%s:%d %.1f != DPMSimple[%u][%u] = %.1f", \
							  __FILE__, __LINE__, x, i, j, g_DPISimple[i][j]); \
							} \
						}

#else

#define cmpm(i, j, x)	/* empty */
#define cmpd(i, j, x)	/* empty */
#define cmpi(i, j, x)	/* empty */

#endif

static void AllocSave(unsigned LA, unsigned LB)
	{
#if	CMP_SIMPLE
	GetSimpleDPMxs(&g_DPMSimpleMx, &g_DPDSimpleMx, &g_DPISimpleMx);
	g_DPMSimple = g_DPMSimpleMx->GetData();
	g_DPDSimple = g_DPDSimpleMx->GetData();
	g_DPISimple = g_DPISimpleMx->GetData();
#endif
	g_MxDPM.Alloc("FastM", LA+1, LB+1);
	g_MxDPD.Alloc("FastD", LA+1, LB+1);
	g_MxDPI.Alloc("FastI", LA+1, LB+1);

	g_MxTBM.Alloc("FastTBM", LA+1, LB+1);
	g_MxTBD.Alloc("FastTBD", LA+1, LB+1);
	g_MxTBI.Alloc("FastTBI", LA+1, LB+1);

	g_DPM = g_MxDPM.GetData();
	g_DPD = g_MxDPD.GetData();
	g_DPI = g_MxDPI.GetData();

	g_TBM = g_MxTBM.GetData();
	g_TBD = g_MxTBD.GetData();
	g_TBI = g_MxTBI.GetData();
	}

static void SAVE_DPM(unsigned i, unsigned j, float x)
	{
	g_DPM[i][j] = x;
#if	CMP_SIMPLE
	if (i > 0 && j > 0)
	asserta(feq(x, g_DPMSimple[i][j]));
#endif
	}

static void SAVE_DPD(unsigned i, unsigned j, float x)
	{
	g_DPD[i][j] = x;
#if	CMP_SIMPLE
	if (i > 0 && j > 0)
	asserta(feq(x, g_DPDSimple[i][j]));
#endif
	}

static void SAVE_DPI(unsigned i, unsigned j, float x)
	{
	g_DPI[i][j] = x;
#if	CMP_SIMPLE
	if (i > 0 && j > 0)
	asserta(feq(x, g_DPISimple[i][j]));
#endif
	}

static void SAVE_TBM(unsigned i, unsigned j, char x)
	{
	g_TBM[i][j] = x;
	}

static void SAVE_TBD(unsigned i, unsigned j, char x)
	{
	g_TBD[i][j] = x;
	}

static void SAVE_TBI(unsigned i, unsigned j, char x)
	{
	g_TBI[i][j] = x;
	}

void GetFastMxs(Mx<float> **M, Mx<float> **D, Mx<float> **I)
	{
	*M = &g_MxDPM;
	*D = &g_MxDPD;
	*I = &g_MxDPI;
	}

#else	// SAVE_FAST

#define	SAVE_DPM(i, j, x)	/* empty */
#define	SAVE_DPD(i, j, x)	/* empty */
#define	SAVE_DPI(i, j, x)	/* empty */

#define	SAVE_TBM(i, j, x)	/* empty */
#define	SAVE_TBD(i, j, x)	/* empty */
#define	SAVE_TBI(i, j, x)	/* empty */

#define AllocSave(LA, LB)	/* empty */

#define cmpm(i, j, x)	/* empty */
#define cmpd(i, j, x)	/* empty */
#define cmpi(i, j, x)	/* empty */

#endif	// SAVE_FAST

float ViterbiFast(const byte *A, unsigned LA, const byte *B, unsigned LB,
  const AlnParams &AP, PathData &PD)
	{
	if (LA*LB > 100*1000*1000)
		Die("ViterbiFast, too long LA=%u, LB=%u", LA, LB);

	AllocBit(LA, LB);
	AllocSave(LA, LB);
	
	StartTimer(ViterbiFast);

	const float * const *Mx = AP.SubstMx;
	float OpenA = AP.LOpenA;
	float ExtA = AP.LExtA;

	byte **TB = g_TBBit;
	float *Mrow = g_DPRow1;
	float *Drow = g_DPRow2;

// Use Mrow[-1], so...
	Mrow[-1] = MINUS_INFINITY;
	for (unsigned j = 0; j <= LB; ++j)
		{
		Mrow[j] = MINUS_INFINITY;
		SAVE_DPM(0, j, MINUS_INFINITY);
		SAVE_TBM(0, j, '?');

		Drow[j] = MINUS_INFINITY;
		SAVE_DPD(0, j, MINUS_INFINITY);
		SAVE_TBD(0, j, '?');
		}
	
// Main loop
	float M0 = float (0);
	SAVE_DPM(0, 0, 0);
	for (unsigned i = 0; i < LA; ++i)
		{
		byte a = A[i];
		const float *MxRow = Mx[a];
		float OpenB = AP.LOpenB;
		float ExtB = AP.LExtB;
		float I0 = MINUS_INFINITY;

		SAVE_TBM(i, 0, '?');

		SAVE_DPI(i, 0, MINUS_INFINITY);
		SAVE_DPI(i, 1, MINUS_INFINITY);

		SAVE_TBI(i, 0, '?');
		SAVE_TBI(i, 1, '?');
		
		byte *TBrow = TB[i];
		for (unsigned j = 0; j < LB; ++j)
			{
			byte b = B[j];
			byte TraceBits = 0;
			float SavedM0 = M0;

		// MATCH
			{
		// M0 = DPM[i][j]
		// I0 = DPI[i][j]
		// Drow[j] = DPD[i][j]
			cmpm(i, j, M0);
			cmpd(i, j, Drow[j]);
			cmpi(i, j, I0);

			float xM = M0;
			SAVE_TBM(i+1, j+1, 'M');
			if (Drow[j] > xM)
				{
				xM = Drow[j];
				TraceBits = TRACEBITS_DM;
				SAVE_TBM(i+1, j+1, 'D');
				}
			if (I0 > xM)
				{
				xM = I0;
				TraceBits = TRACEBITS_IM;
				SAVE_TBM(i+1, j+1, 'I');
				}
			M0 = Mrow[j];
			cmpm(i, j+1, M0);

			Mrow[j] = xM + MxRow[b];
		// Mrow[j] = DPM[i+1][j+1])
			SAVE_DPM(i+1, j+1, Mrow[j]);
			}
			
		// DELETE
			{
		// SavedM0 = DPM[i][j]
		// Drow[j] = DPD[i][j]
			cmpm(i, j, SavedM0);
			cmpd(i, j, Drow[j]);

			float md = SavedM0 + OpenB;
			Drow[j] += ExtB;
			SAVE_TBD(i+1, j, 'D');
			if (md >= Drow[j])
				{
				Drow[j] = md;
				TraceBits |= TRACEBITS_MD;
				SAVE_TBD(i+1, j, 'M');
				}
		// Drow[j] = DPD[i+1][j]
			SAVE_DPD(i+1, j, Drow[j]);
			}
			
		// INSERT
			{
		// SavedM0 = DPM[i][j]
		// I0 = DPI[i][j]
			cmpm(i, j, SavedM0);
			cmpi(i, j, I0);
			
			float mi = SavedM0 + OpenA;
			I0 += ExtA;
			SAVE_TBI(i, j+1, 'I');
			if (mi >= I0)
				{
				I0 = mi;
				TraceBits |= TRACEBITS_MI;
				SAVE_TBI(i, j+1, 'M');
				}
		// I0 = DPI[i][j+1]
			SAVE_DPI(i, j+1, I0);
			}
			
			OpenB = AP.OpenB;
			ExtB = AP.ExtB;
			
			TBrow[j] = TraceBits;
			}
		
	// Special case for end of Drow[]
		{
	// M0 = DPM[i][LB]
	// Drow[LB] = DPD[i][LB]
		
		TBrow[LB] = 0;
		float md = M0 + AP.ROpenB;
		Drow[LB] += AP.RExtB;
		SAVE_TBD(i+1, LB, 'D');
		if (md >= Drow[LB])
			{
			Drow[LB] = md;
			TBrow[LB] = TRACEBITS_MD;
			SAVE_TBD(i+1, LB, 'M');
			}
	// Drow[LB] = DPD[i+1][LB]
		SAVE_DPD(i+1, LB, Drow[LB]);
		}
		
		SAVE_DPM(i+1, 0, MINUS_INFINITY);
		M0 = MINUS_INFINITY;

		OpenA = AP.OpenA;
		ExtA = AP.ExtA;
		}
	
	SAVE_TBM(LA, 0, '?');

// Special case for last row of DPI
	byte *TBrow = TB[LA];
	float I1 = MINUS_INFINITY;

	SAVE_DPI(LA, 0, MINUS_INFINITY);
	SAVE_TBI(LA, 0, '?');

	SAVE_DPI(LA, 1, MINUS_INFINITY);
	SAVE_TBI(LA, 1, '?');

	for (unsigned j = 1; j < LB; ++j)
		{
	// Mrow[j-1] = DPM[LA][j]
	// I1 = DPI[LA][j]
		
		TBrow[j] = 0;
		float mi = Mrow[int(j)-1] + AP.ROpenA;
		I1 += AP.RExtA;
		SAVE_TBI(LA, j+1, 'I');
		if (mi > I1)
			{
			I1 = mi;
			TBrow[j] = TRACEBITS_MI;
			SAVE_TBI(LA, j+1, 'M');
			}
		SAVE_DPI(LA, j+1, I1);
		}
	
	float FinalM = Mrow[LB-1];
	float FinalD = Drow[LB];
	float FinalI = I1;
// FinalM = DPM[LA][LB]
// FinalD = DPD[LA][LB]
// FinalI = DPI[LA][LB]
	
	float Score = FinalM;
	byte State = 'M';
	if (FinalD > Score)
		{
		Score = FinalD;
		State = 'D';
		}
	if (FinalI > Score)
		{
		Score = FinalI;
		State = 'I';
		}

	EndTimer(ViterbiFast);
	TraceBackBit(LA, LB, State, PD);

#if	SAVE_FAST
	g_MxDPM.LogMe();
	g_MxDPD.LogMe();
	g_MxDPI.LogMe();

	g_MxTBM.LogMe();
	g_MxTBD.LogMe();
	g_MxTBI.LogMe();
#endif

	return Score;
	}
