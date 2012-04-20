#if	UCHIMES

#include "myutils.h"
#include "chime.h"
#include "ultra.h"
#include <set>

const float MAX_WORD_COUNT_DROP = 1;

void SortDescending(const vector<float> &Values, vector<unsigned> &Order);
bool GlobalAlign(const SeqData &Query, const SeqData &Target, string &Path);
double GetFractIdGivenPath(const byte *A, const byte *B, const char *Path);
void USort(const SeqData &Query, const SeqDB &DB, vector<float> &WordCounts,
  vector<unsigned> &Order);

void AddTargets(SeqDB &DB, const SeqData &Query, set<unsigned> &TargetIndexes)
	{
	const unsigned SeqCount = DB.GetSeqCount();
	if (SeqCount == 0)
		return;

	vector<float> WordCounts;
	vector<unsigned> Order;
	USort(Query, DB, WordCounts, Order);
	asserta(SIZE(Order) == SeqCount);
	unsigned TopSeqIndex = Order[0];
	float TopWordCount = WordCounts[TopSeqIndex];
	for (unsigned i = 0; i < SeqCount; ++i)
		{
		unsigned SeqIndex = Order[i];
		float WordCount = WordCounts[SeqIndex];
		if (TopWordCount - WordCount > MAX_WORD_COUNT_DROP)
			return;
		TargetIndexes.insert(SeqIndex);
		}
	}

#endif
