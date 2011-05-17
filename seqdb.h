//uchime by Robert C. Edgar http://drive5.com/uchime This code is donated to the public domain.

#ifndef seqdb_h
#define seqdb_h

#include <vector>
#include <map>
#include "myutils.h"

struct SeqData;

using namespace std;

struct SeqDB
	{
private:
	SeqDB(const SeqDB &rhs);
	SeqDB &operator=(const SeqDB &rhs);

public:
	string m_FileName;
	char **m_Labels;
	byte **m_Seqs;
	unsigned *m_SeqLengths;
	unsigned m_SeqCount;
	unsigned m_Size;

	bool m_Aligned;
	bool m_IsNucleo;
	bool m_IsNucleoSet;

public:
	SeqDB();
	~SeqDB();
	void Clear(bool ctor = false);
	void InitEmpty(bool Nucleo);

	unsigned AddSeq(const char *Label, const byte *Seq, unsigned L);

	byte *GetSeq(unsigned SeqIndex) const
		{
		asserta(SeqIndex < m_SeqCount);
		return m_Seqs[SeqIndex];
		}

	const char *GetLabel(unsigned SeqIndex) const
		{
		asserta(SeqIndex < m_SeqCount);
		return m_Labels[SeqIndex];
		}

	unsigned GetSeqLength(unsigned SeqIndex) const
		{
		asserta(SeqIndex < m_SeqCount);
		return m_SeqLengths[SeqIndex];
		}

	unsigned GetSeqCount() const
		{
		return m_SeqCount;
		}

	unsigned GetPairCount() const
		{
		unsigned SeqCount = GetSeqCount();
		return (SeqCount*(SeqCount - 1))/2;
		}

	unsigned GetPairIndex(unsigned SeqIndex1, unsigned SeqIndex2) const
		{
		if (SeqIndex1 > SeqIndex2)
			return (SeqIndex1*(SeqIndex1 - 1))/2 + SeqIndex2;
		return (SeqIndex2*(SeqIndex2 - 1))/2 + SeqIndex1;
		}

	unsigned GetColCount() const
		{
		if (!m_Aligned)
			Die("SeqDB::GetColCount, not aligned");
		if (m_SeqCount == 0)
			Die("SeqDB::GetColCount, empty");
		return m_SeqLengths[0];
		}

	bool IsNucleo() const
		{
		asserta(m_IsNucleoSet);
		return m_IsNucleo;
		}

	void GetSeqData(unsigned Id, SeqData &Buffer) const;

	unsigned GetMaxLabelLength() const;
	unsigned GetMaxSeqLength() const;
	void SetIsNucleo();
	unsigned GetIndex(const char *Label) const;
	void MakeLabelToIndex(map<string, unsigned> &LabelToIndex);

	void LogMe() const;
	void FromFasta(const string &FileName, bool AllowGaps = false);

	void ToFasta(const string &FileName) const;
	void ToFasta(FILE *f, unsigned SeqIndex) const;
	void SeqToFasta(FILE *f, unsigned SeqIndex, bool WithLabel = false) const;

	unsigned GetTotalLength() const;
	};

bool isgap(byte c);

#endif
