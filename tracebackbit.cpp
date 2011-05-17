//uchime by Robert C. Edgar http://drive5.com/uchime This code is donated to the public domain.

#include "dp.h"

#define TRACE	0

Mx<byte> g_Mx_TBBit;
byte **g_TBBit;
float *g_DPRow1;
float *g_DPRow2;
static float *g_DPBuffer1;
static float *g_DPBuffer2;

static unsigned g_CacheLB;

void AllocBit(unsigned LA, unsigned LB)
	{
	g_Mx_TBBit.Alloc("TBBit", LA+1, LB+1);
	g_TBBit = g_Mx_TBBit.GetData();
	if (LB > g_CacheLB)
		{
		MYFREE(g_DPBuffer1, g_CacheLB, AllocBit);
		MYFREE(g_DPBuffer2, g_CacheLB, AllocBit);

		g_CacheLB = LB + 128;

	// Allow use of [-1]
		//g_DPBuffer1 = myalloc<float>(g_CacheLB+3);
		//g_DPBuffer2 = myalloc<float>(g_CacheLB+3);
		g_DPBuffer1 = MYALLOC(float, g_CacheLB+3, AllocBit);
		g_DPBuffer2 = MYALLOC(float, g_CacheLB+3, AllocBit);
		g_DPRow1 = g_DPBuffer1 + 1;
		g_DPRow2 = g_DPBuffer2 + 1;
		}
	}

void TraceBackBit(unsigned LA, unsigned LB, char State, PathData &PD)
	{
	PD.Alloc(LA+LB);

	StartTimer(TraceBackBit);
	char *PathPtr = PD.Back;
	*PathPtr = 0;

	byte **TB = g_TBBit;

#if	TRACE
	Log("\n");
	Log("TraceBackBit\n");
#endif

	size_t i = LA;
	size_t j = LB;
	for (;;)
		{
#if	TRACE
		Log("i=%3d  j=%3d  state=%c\n", (int) i, (int) j, State);
#endif
		if (i == 0 && j == 0)
			break;

		--PathPtr;
		*PathPtr = State;

		byte t;
		switch (State)
			{
		case 'M':
			asserta(i > 0 && j > 0);
			t = TB[i-1][j-1];
			if (t & TRACEBITS_DM)
				State = 'D';
			else if (t & TRACEBITS_IM)
				State = 'I';
			else
				State = 'M';
			--i;
			--j;
			break;
		case 'D':
			asserta(i > 0);
			t = TB[i-1][j];
			if (t & TRACEBITS_MD)
				State = 'M';
			else
				State = 'D';
			--i;
			break;

		case 'I':
			asserta(j > 0);
			t = TB[i][j-1];
			if (t & TRACEBITS_MI)
				State = 'M';
			else
				State = 'I';
			--j;
			break;

		default:
			Die("TraceBackBit, invalid state %c", State);
			}
		}
	PD.Start = PathPtr;
	EndTimer(TraceBackBit);
	}

void TraceBackBitSW(unsigned LA, unsigned LB, unsigned Besti, unsigned Bestj,
  unsigned &Leni, unsigned &Lenj, PathData &PD)
	{
	PD.Alloc(LA+LB);

	StartTimer(TraceBackBitSW);
	char *PathPtr = PD.Back;
	*PathPtr = 0;

	byte **TB = g_TBBit;

#if	TRACE
	Log("\n");
	Log("TraceBackBitSW\n");
#endif

	unsigned i = Besti;
	unsigned j = Bestj;
	char State = 'M';
	for (;;)
		{
#if	TRACE
		Log("i=%3d  j=%3d  state=%c\n", (int) i, (int) j, State);
#endif
		--PathPtr;
		*PathPtr = State;

		byte t;
		switch (State)
			{
		case 'M':
			asserta(i > 0 && j > 0);
			t = TB[i-1][j-1];
			if (t & TRACEBITS_DM)
				State = 'D';
			else if (t & TRACEBITS_IM)
				State = 'I';
			else if (t & TRACEBITS_SM)
				{
				Leni = Besti - i + 1;
				Lenj = Bestj - j + 1;
				PD.Start = PathPtr;
				EndTimer(TraceBackBitSW);
				return;
				}
			else
				State = 'M';
			--i;
			--j;
			break;
		case 'D':
			asserta(i > 0);
			t = TB[i-1][j];
			if (t & TRACEBITS_MD)
				State = 'M';
			else
				State = 'D';
			--i;
			break;

		case 'I':
			asserta(j > 0);
			t = TB[i][j-1];
			if (t & TRACEBITS_MI)
				State = 'M';
			else
				State = 'I';
			--j;
			break;

		default:
			Die("TraceBackBitSW, invalid state %c", State);
			}
		}
	}
