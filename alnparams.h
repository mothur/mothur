#ifndef alnparams_h
#define alnparams_h

struct HSPData;

// Gap penalty scores are negative
// (i.e., are scores, not penalties).
struct AlnParams
	{
	const char *SubstMxName;
	const float * const *SubstMx;

	bool Nucleo;
	bool NucleoSet;

// Local gaps
	float LocalOpen;
	float LocalExt;

// Global internal gaps
	float OpenA;
	float OpenB;

	float ExtA;
	float ExtB;

// Global terminal gaps
	float LOpenA;
	float LOpenB;
	float ROpenA;
	float ROpenB;

	float LExtA;
	float LExtB;
	float RExtA;
	float RExtB;

	void Clear();
	void SetLocal(float Open, float Ext);
	void Init2(const float * const *Mx, float Open, float Ext);
	void Init4(const float * const *Mx, float Open, float Ext, float TermOpen, float TermExt);
	void Init(const AlnParams &AP, const HSPData &HSP, unsigned LA, unsigned LB);
	void InitFromCmdLine(bool Nucleo);
	void SetMxFromCmdLine(bool Nucleo);
	void SetPenalties(const string &OpenStr, const string &ExtStr);
	float GetLocalOpen() const;
	float GetLocalExt() const;
	bool GetIsNucleo() const;

	bool Is2() const;
	bool Is4() const;
	const char *GetType() const;

	void LogMe() const;
	};

const float OBVIOUSLY_WRONG_PENALTY = 1000.0;

#endif // alnparams_h
