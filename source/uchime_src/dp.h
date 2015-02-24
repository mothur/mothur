#ifndef dp_h
#define dp_h

#define SAVE_FAST	0

#include "myutils.h"
#include "mx.h"
#include "seqdb.h"
#include "diagbox.h"
#include "path.h"
#include "alnparams.h"
#include "alnheuristics.h"
#include "hspfinder.h"

typedef void (*OnPathFn)(const string &Path, bool Full);

enum XType
	{
	XType_Full=1,
	XType_Fwd=2,
	XType_Bwd=3,
	};

// public
float ViterbiBrute(const byte *A, unsigned LA, const byte *B, unsigned LB, 
  unsigned DiagLo, unsigned DiagHi, const AlnParams &AP, PathData &PD);

float ViterbiSimple(const byte *A, unsigned LA, const byte *B, unsigned LB,
  const AlnParams &AP, PathData &PD);

float ViterbiSimpleBand(const byte *A, unsigned LA, const byte *B, unsigned LB,
  const AlnParams &AP, unsigned DiagLo, unsigned DiagHi, PathData &PD);

float ViterbiFast(const byte *A, unsigned LA, const byte *B, unsigned LB,
  const AlnParams &AP, PathData &PD);

float ViterbiFastBand(const byte *A, unsigned LA, const byte *B, unsigned LB,
  unsigned DiagLo, unsigned DiagHi, const AlnParams &AP, PathData &PD);

float ViterbiFastMainDiag(const byte *A, unsigned LA, const byte *B, unsigned LB,
  unsigned BandRadius, const AlnParams &AP, PathData &PD);

float XDropFwdSimple(const byte *A, unsigned LA, const byte *B, unsigned LB,
  const AlnParams &AP, float XDrop, unsigned &Leni, unsigned &Lenj, PathData &PD);

float XDropBwdSimple(const byte *A, unsigned LA, const byte *B, unsigned LB,
  const AlnParams &AP, float XDrop, unsigned &Leni, unsigned &Lenj, PathData &PD);

float XDropFwdFast(const byte *A, unsigned LA, const byte *B, unsigned LB,
  const AlnParams &AP, float XDrop, unsigned &Leni, unsigned &Lenj, PathData &PD);

float XDropBwdFast(const byte *A, unsigned LA, const byte *B, unsigned LB,
  const AlnParams &AP, float XDrop, unsigned &Leni, unsigned &Lenj, PathData &PD);

void XDropAlign(const byte *A, unsigned LA, const byte *B, unsigned LB,
  unsigned AncLoi, unsigned AncLoj, unsigned AncLen, const AlnParams &AP,
  float XDrop, HSPData &HSP, PathData &PD);

float SWSimple(const byte *A, unsigned LA, const byte *B, unsigned LB,
  const AlnParams &AP, unsigned &Loi, unsigned &Leni, unsigned &Lenj,
  unsigned &Hij, PathData &PD);

float SWFast(const byte *A, unsigned LA, const byte *B, unsigned LB,
  const AlnParams &AP, unsigned &Loi, unsigned &Leni, unsigned &Lenj,
  unsigned &Hij, PathData &PD);

void SWFast2(const SeqData &SA, const SeqData &SB, const AlnParams &AP,
  HSPData &HSP, PathData &PD);

void SWSimple2(const SeqData &SA, const SeqData &SB, const AlnParams &AP,
  HSPData &HSP, PathData &PD);

float SWUngapped(const byte *A, unsigned LA, const byte *B, unsigned LB,
  const float * const *SubstMx, unsigned &LoA, unsigned &LoB, unsigned &Len);

void SWUngapped2(const SeqData &SA, const SeqData &SB, const AlnParams &AP,
  HSPData &HSP);

float SWFastNTB(const byte *A, unsigned LA, const byte *B, unsigned LB,
  const AlnParams &AP);

void GlobalAlignBand(const byte *A, unsigned LA, const byte *B, unsigned LB,
  const AlnParams &AP, unsigned BandRadius, PathData &PD);

bool GlobalAlign(const SeqData &Query, const SeqData &Target, const AlnParams &AP,
  const AlnHeuristics &AH, HSPFinder &HF, float MinFractId, float &HSPFractId,
  PathData &PD);

bool GlobalAlign(const SeqData &Query, const SeqData &Target, string &Path);

void GetBruteMxs(Mx<float> **M, Mx<float> **D, Mx<float> **I);
void GetSimpleDPMxs(Mx<float> **M, Mx<float> **D, Mx<float> **I);
void GetSimpleBandMxs(Mx<float> **M, Mx<float> **D, Mx<float> **I);
void GetXDropFwdSimpleDPMxs(Mx<float> **M, Mx<float> **D, Mx<float> **I);
#if	SAVE_FAST
void GetFastMxs(Mx<float> **M, Mx<float> **D, Mx<float> **I);
void GetFastBandMxs(Mx<float> **M, Mx<float> **D, Mx<float> **I);
#endif

// private
void TraceBackBit(unsigned LA, unsigned LB, char State, PathData &PD);
void TraceBackBitSW(unsigned LA, unsigned LB, unsigned Besti, unsigned Bestj,
  unsigned &Leni, unsigned &Lenj, PathData &PD);
void EnumPaths(unsigned L1, unsigned L2, bool SubPaths, OnPathFn OnPath);
void AllocBit(unsigned LA, unsigned LB);

const byte TRACEBITS_DM = 0x01;
const byte TRACEBITS_IM = 0x02;
const byte TRACEBITS_MD = 0x04;
const byte TRACEBITS_MI = 0x08;
const byte TRACEBITS_SM = 0x10;
const byte TRACEBITS_UNINIT = ~0x1f;

extern Mx<byte> g_Mx_TBBit;
extern float *g_DPRow1;
extern float *g_DPRow2;
extern byte **g_TBBit;

static inline void Max_xM(float &Score, float MM, float DM, float IM, byte &State)
	{
	Score = MM;
	State = 'M';

	if (DM > Score)
		{
		Score = DM;
		State = 'D';
		}
	if (IM > Score)
		{
		Score = IM;
		State = 'I';
		}
	}

static inline void Max_xD(float &Score, float MD, float DD, byte &State)
	{
	if (MD >= DD)
		{
		Score = MD;
		State = 'M';
		}
	else
		{
		Score = DD;
		State = 'D';
		}
	}

static inline void Max_xI(float &Score, float MI, float II, byte &State)
	{
	if (MI >= II)
		{
		Score = MI;
		State = 'M';
		}
	else
		{
		Score = II;
		State = 'I';
		}
	}

#endif // dp_h
