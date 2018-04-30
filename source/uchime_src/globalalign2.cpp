#if	UCHIMES

#include "dp.h"
#include "seq.h"

static AlnParams g_AP;
static bool g_APInitDone = false;

bool GlobalAlign(const SeqData &Query, const SeqData &Target, PathData &PD)
	{
	if (!g_APInitDone)
		{
		g_AP.InitFromCmdLine(true);
		g_APInitDone = true;
		}

	ViterbiFast(Query.Seq, Query.L, Target.Seq, Target.L, g_AP, PD);
	return true;
	}

bool GlobalAlign(const SeqData &Query, const SeqData &Target, string &Path)
	{
	PathData PD;
	GlobalAlign(Query, Target, PD);
	Path = string(PD.Start);
	return true;
	}

bool GlobalAlign(const SeqData &Query, const SeqData &Target, const AlnParams &/*AP*/,
  const AlnHeuristics &AH, HSPFinder &/*HF*/, float /*MinFractId*/, float &/*HSPId*/, PathData &PD)
	{
	PD.Clear();
	string Path;
	bool Found = GlobalAlign(Query, Target, Path);
	if (!Found)
		return false;
	unsigned n = SIZE(Path);
	PD.Alloc(n+1);
	memcpy(PD.Front, Path.c_str(), n);
	PD.Start = PD.Front;
	PD.Start[n] = 0;
	return true;
	}

#endif // UCHIMES
