#include "myutils.h"
#include "chime.h"
#include "ultra.h"
#include <set>

void AddTargets(Ultra &U, const SeqData &Query, set<unsigned> &TargetIndexes);

void GetChunkInfo(unsigned L, unsigned &Length, vector<unsigned> &Los)
	{
	Los.clear();

	if (L <= opt_minchunk)
		{
		Length = L;
		Los.push_back(0);
		return;
		}

	Length = (L - 1)/opt_chunks + 1;
	if (Length < opt_minchunk)
		Length = opt_minchunk;

	unsigned Lo = 0;
	for (;;)
		{
		if (Lo + Length >= L)
			{
			Lo = L - Length - 1;
			Los.push_back(Lo);
			return;
			}
		Los.push_back(Lo);
		Lo += Length;
		}
	}

void GetCandidateParents(Ultra &U, const SeqData &QSD, float AbQ,
  vector<unsigned> &Parents)
	{
	Parents.clear();

	set<unsigned> TargetIndexes;

	unsigned QL = QSD.L;

	SeqData QuerySD = QSD;

	unsigned ChunkLength;
	vector<unsigned> ChunkLos;
	GetChunkInfo(QL, ChunkLength, ChunkLos);
	unsigned ChunkCount = SIZE(ChunkLos);
	for (unsigned ChunkIndex = 0; ChunkIndex < ChunkCount; ++ChunkIndex)
		{
		unsigned Lo = ChunkLos[ChunkIndex];
		asserta(Lo + ChunkLength <= QL);

		const byte *Chunk = QSD.Seq + Lo;

	// THIS MESSES UP --self!!
		//char Prefix[32];
		//sprintf(Prefix, "%u|", Lo);
		//string ChunkLabel = string(Prefix) + string(QSD.Label);

		//QuerySD.Label = ChunkLabel.c_str();
		QuerySD.Seq = Chunk;
		QuerySD.L = ChunkLength;

		AddTargets(U, QuerySD, TargetIndexes);

		Lo += ChunkLength;
		}

	for (set<unsigned>::const_iterator p = TargetIndexes.begin();
	  p != TargetIndexes.end(); ++p)
		{
		unsigned TargetIndex = *p;
		bool Accept = true;
		if (AbQ > 0.0f)
			{
			const char *TargetLabel = U.GetSeedLabel(TargetIndex);
			float AbT = GetAbFromLabel(string(TargetLabel));
			if (AbT > 0.0f && AbT < opt_abskew*AbQ)
				Accept = false;
			}

		if (Accept)
			Parents.push_back(TargetIndex);
		}
	}
