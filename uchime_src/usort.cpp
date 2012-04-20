#if	UCHIMES

#include "myutils.h"
#include "seqdb.h"
#include "seq.h"
#include "alpha.h"

void SortDescending(const vector<float> &Values, vector<unsigned> &Order);

static byte *g_QueryHasWord;
static unsigned g_WordCount;

unsigned GetWord(const byte *Seq)
	{
	unsigned Word = 0;
	const byte *Front = Seq;
	for (unsigned i = 0; i < opt_w; ++i)
		{
		unsigned Letter = g_CharToLetterNucleo[*Front++];
		Word = (Word*4) + Letter;
		}
	return Word;
	}

static void SetQuery(const SeqData &Query)
	{
	if (g_QueryHasWord == 0)
		{
		g_WordCount = 4;
		for (unsigned i = 1; i < opt_w; ++i)
			g_WordCount *= 4;

		g_QueryHasWord = myalloc(byte, g_WordCount);
		}

	memset(g_QueryHasWord, 0, g_WordCount);

	if (Query.L <= opt_w)
		return;

	const unsigned L = Query.L - opt_w + 1;
	const byte *Seq = Query.Seq;
	for (unsigned i = 0; i < L; ++i)
		{
		unsigned Word = GetWord(Seq++);
		g_QueryHasWord[Word] = 1;
		}
	}

static unsigned GetUniqueWordsInCommon(const SeqData &Target)
	{
	if (Target.L <= opt_w)
		return 0;

	unsigned Count = 0;
	const unsigned L = Target.L - opt_w + 1;
	const byte *Seq = Target.Seq;
	for (unsigned i = 0; i < L; ++i)
		{
		unsigned Word = GetWord(Seq++);
		if (g_QueryHasWord[Word])
			++Count;
		}
	return Count;
	}

void USort(const SeqData &Query, const SeqDB &DB, vector<float> &WordCounts, 
  vector<unsigned> &Order)
	{
	WordCounts.clear();
	Order.clear();

	SetQuery(Query);

	const unsigned SeqCount = DB.GetSeqCount();
	for (unsigned SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		SeqData Target;
		DB.GetSeqData(SeqIndex, Target);
		float WordCount = (float) GetUniqueWordsInCommon(Target);
		WordCounts.push_back(WordCount);
		}
	SortDescending(WordCounts, Order);
	}

#endif // UCHIMES
