#ifndef chime_h
#define chime_h

#include "seq.h"

struct ChimeHit2
	{
	string QLabel;
	string ALabel;
	string BLabel;
	string Q3;
	string A3;
	string B3;

	//unsigned LY, LN, LA, LD;
	//unsigned RY, RN, RA, RD;
	double PctIdQT, PctIdQA, PctIdQB, PctIdQM, PctIdAB;

	unsigned ColLo;
	unsigned ColXLo;
	unsigned ColXHi;
	unsigned ColHi;
	unsigned QXLo;
	unsigned QXHi;

	double Div;
	double Score;
	double H;

	unsigned CS_LY, CS_LN, CS_LA, CS_RY, CS_RN, CS_RA;

	float AbQ;
	float AbA;
	float AbB;

	ChimeHit2()
		{
		Clear();
		}

	void Clear()
		{
		Q3.clear();
		A3.clear();
		B3.clear();
		QLabel.clear();
		ALabel.clear();
		BLabel.clear();

		//LY = LN = LA = LD = UINT_MAX;
		//RY = RN = RA = RD = UINT_MAX;
		ColLo = ColHi = QXLo = QXHi = ColXLo = ColXHi = UINT_MAX;
		CS_LY = CS_LN = CS_LA = CS_RY = CS_RN = CS_RA = UINT_MAX;
		PctIdQT = PctIdQA = PctIdQB = PctIdQM = PctIdAB = -1.0;
		Div = -1.0;
		H = -1.0;
		Score = -1.0;
		AbQ = AbA = AbB = -1.0f;
		};

	bool Accept() const
		{
		return Score >= opt_minh && Div >= opt_mindiv && CS_LY >= opt_mindiffs && CS_RY >= opt_mindiffs;
		}

	void LogMe() const
		{
		Log("@L %c ", yon(Score >= 1.0 && Div >= 1.0));
		Log(" %.4f", Score);
		Log(" LY %u LN %u LA %u", CS_LY, CS_LN, CS_LA);
		Log(" RY %u RN %u RA %u", CS_RY, CS_RN, CS_RA);
		Log(" Div %.1f%%", Div);
		Log(" Q=%s", QLabel.c_str());
		Log(" A=%s", ALabel.c_str());
		Log(" B=%s", BLabel.c_str());
		Log(" QA %.1f%% QB=%.1f%% AB=%.1f%% QM=%.1f%%", PctIdQA, PctIdQB, PctIdAB, PctIdQM);
		Log("\n");
		}

	bool operator<(const ChimeHit2 &rhs) const
		{
		if (Score == rhs.Score)
			return Div > rhs.Div;
		return Score > rhs.Score;
		}
	};

static inline bool isacgt(char c)
	{
	return c == 'A' || c == 'C' || c == 'G' || c == 'T';
	}

static bool inline isgap(char c)
	{
	return c == '-' || c == '.';
	}

void GetChunkInfo(unsigned L, unsigned &Length, vector<unsigned> &Los);
float GetAbFromLabel(const string &Label);
void WriteChimeHitCS(FILE *f, const ChimeHit2 &Hit);
void WriteChimeHit(FILE *f, const ChimeHit2 &Hit);
void WriteChimeFileHdr(FILE *f);

#endif // chime_h
