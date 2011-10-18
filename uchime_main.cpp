#include "myutils.h"
#include "chime.h"
#include "seqdb.h"
#include "dp.h"
#include "ultra.h"
#include "hspfinder.h"
#include <algorithm>
#include <set>

bool SearchChime(Ultra &U, const SeqData &QSD, float QAb, 
  const AlnParams &AP, const AlnHeuristics &AH, HSPFinder &HF,
  float MinFractId, ChimeHit2 &Hit);

FILE *g_fUChime;
FILE *g_fUChimeAlns;
const vector<float> *g_SortVecFloat;
bool g_UchimeDeNovo = false;

void Usage()
	{
	printf("\n");
	printf("UCHIME %s by Robert C. Edgar\n", MY_VERSION);
	printf("http://www.drive5.com/uchime\n");
	printf("\n");
	printf("This software is donated to the public domain\n");
	printf("\n");

	printf(
#include "help.h"
		);
	}

void SetBLOSUM62()
	{
	Die("SetBLOSUM62 not implemented");
	}

void ReadSubstMx(const string &/*FileName*/, Mx<float> &/*Mxf*/)
	{
	Die("ReadSubstMx not implemented");
	}

void LogAllocs()
	{
	/*empty*/
	}

static bool CmpDescVecFloat(unsigned i, unsigned j)
	{
	return (*g_SortVecFloat)[i] > (*g_SortVecFloat)[j];
	}

void Range(vector<unsigned> &v, unsigned N)
	{
	v.clear();
	v.reserve(N);
	for (unsigned i = 0; i < N; ++i)
		v.push_back(i);
	}

void SortDescending(const vector<float> &Values, vector<unsigned> &Order)
	{
	StartTimer(Sort);
	const unsigned N = SIZE(Values);
	Range(Order, N);
	g_SortVecFloat = &Values;
	sort(Order.begin(), Order.end(), CmpDescVecFloat);
	EndTimer(Sort);
	}

float GetAbFromLabel(const string &Label)
	{
	vector<string> Fields;
	Split(Label, Fields, '/');
	const unsigned N = SIZE(Fields);
	for (unsigned i = 0; i < N; ++i)
		{
		const string &Field = Fields[i];
		if (Field.substr(0, 3) == "ab=")
			{
			string a = Field.substr(3, string::npos);
			return (float) atof(a.c_str());
			}
		}
	if (g_UchimeDeNovo)
		Die("Missing abundance /ab=xx/ in label >%s", Label.c_str());
	return 0.0;
	}

int uchime_main(int argc, char *argv[])
	{
		
	MyCmdLine(argc, argv);

	if (argc < 2)
		{
		Usage();
		return 0;
		}

	if (opt_version)
		{
		printf("uchime v" MY_VERSION ".%s\n", SVN_VERSION);
		return 0;
		}

	printf("uchime v" MY_VERSION ".%s\n", SVN_VERSION);
	printf("by Robert C. Edgar\n");
	printf("http://drive5.com/uchime\n");
	printf("This code is donated to the public domain.\n");
	printf("\n");
	if (!optset_w)
		opt_w = 8;
	
	float MinFractId = 0.95f;
	if (optset_id)
		MinFractId = (float) opt_id;

	Log("%8.2f  minh\n", opt_minh);
	Log("%8.2f  xn\n", opt_xn);
	Log("%8.2f  dn\n", opt_dn);
	Log("%8.2f  xa\n", opt_xa);
	Log("%8.2f  mindiv\n", opt_mindiv);
	Log("%8u  maxp\n", opt_maxp);

	if (opt_input == "" && opt_uchime != "")
		opt_input = opt_uchime;

	if (opt_input == "")
		Die("Missing --input");

	g_UchimeDeNovo = (opt_db == "");

	if (opt_uchimeout != "")
		g_fUChime = CreateStdioFile(opt_uchimeout);

	if (opt_uchimealns != "")
		g_fUChimeAlns = CreateStdioFile(opt_uchimealns);

	SeqDB Input;
	SeqDB DB;

	Input.FromFasta(opt_input);
	if (!Input.IsNucleo())
		Die("Input contains amino acid sequences");

	const unsigned QuerySeqCount = Input.GetSeqCount();
	vector<unsigned> Order;
	for (unsigned i = 0; i < QuerySeqCount; ++i)
		Order.push_back(i);

	if (g_UchimeDeNovo)
		{
		vector<float> Abs;
		for (unsigned i = 0; i < QuerySeqCount; ++i)
			{
			const char *Label = Input.GetLabel(i);
			float Ab = GetAbFromLabel(Label);
			Abs.push_back(Ab);
			}
		SortDescending(Abs, Order);
		DB.m_IsNucleoSet = true;
		DB.m_IsNucleo = true;
		}
	else
		{
		DB.FromFasta(opt_db);
		if (!DB.IsNucleo())
			Die("Database contains amino acid sequences");
		}

	vector<ChimeHit2> Hits;
	unsigned HitCount = 0;
	for (unsigned i = 0; i < QuerySeqCount; ++i)
		{
		unsigned QuerySeqIndex = Order[i];

		SeqData QSD;
		Input.GetSeqData(QuerySeqIndex, QSD);

		float QAb = -1.0;
		if (g_UchimeDeNovo)
			QAb = GetAbFromLabel(QSD.Label);

		ChimeHit2 Hit;
		AlnParams &AP = *(AlnParams *) 0;
		AlnHeuristics &AH = *(AlnHeuristics *) 0;
		HSPFinder &HF = *(HSPFinder *) 0;
		bool Found = SearchChime(DB, QSD, QAb, AP, AH, HF, MinFractId, Hit);
		if (Found)
			++HitCount;
		else
			{
			if (g_UchimeDeNovo)
				DB.AddSeq(QSD.Label, QSD.Seq, QSD.L);
			}

		WriteChimeHit(g_fUChime, Hit);

		ProgressStep(i, QuerySeqCount, "%u/%u chimeras found (%.1f%%)", HitCount, i, Pct(HitCount, i+1));
		}

	Log("\n");
	Log("%s: %u/%u chimeras found (%.1f%%)\n",
	  opt_input.c_str(), HitCount, QuerySeqCount, Pct(HitCount, QuerySeqCount));

	CloseStdioFile(g_fUChime);
	CloseStdioFile(g_fUChimeAlns);

	ProgressExit();
	return 0;
	}
