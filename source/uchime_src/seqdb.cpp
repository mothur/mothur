#include "myutils.h"
#include "seqdb.h"
#include "alpha.h"
#include "timing.h"
#include "sfasta.h"
#include "seq.h"

void SeqToFasta(FILE *f, const char *Label, const byte *Seq, unsigned L)
	{
	const unsigned ROWLEN = 80;
	if (Label != 0)
		fprintf(f, ">%s\n", Label);
	unsigned BlockCount = (L + ROWLEN - 1)/ROWLEN;
	for (unsigned BlockIndex = 0; BlockIndex < BlockCount; ++BlockIndex)
		{
		unsigned From = BlockIndex*ROWLEN;
		unsigned To = From + ROWLEN;
		if (To >= L)
			To = L;
		for (unsigned Pos = From; Pos < To; ++Pos)
			fputc(Seq[Pos], f);
		fputc('\n', f);
		}
	}

SeqDB::~SeqDB()
	{
	Clear();
	}

SeqDB::SeqDB()
	{
	Clear(true);
	}

void SeqDB::Clear(bool ctor)
	{
	if (!ctor)
		{
		for (unsigned i = 0; i < m_SeqCount; ++i)
			{
			unsigned n = strlen(m_Labels[i]);
			MYFREE(m_Labels[i], n, SeqDB);
			MYFREE(m_Seqs[i], m_SeqLengths[i], SeqDB);
			}
		MYFREE(m_Labels, m_Size, SeqDB);
		MYFREE(m_Seqs, m_Size, SeqDB);
		MYFREE(m_SeqLengths, m_Size, SeqDB);
		}

	m_FileName.clear();
	m_SeqCount = 0;
	m_Size = 0;

	m_Labels = 0;
	m_Seqs = 0;
	m_SeqLengths = 0;

	m_Aligned = false;
	m_IsNucleo = false;
	m_IsNucleoSet = false;
	}

void SeqDB::InitEmpty(bool Nucleo)
	{
	Clear();
	m_IsNucleo = Nucleo;
	m_IsNucleoSet = true;
	}

void SeqDB::FromFasta(const string &FileName, bool AllowGaps)
	{
	Clear();
	m_FileName = FileName;
	SFasta SF;

	SF.Open(FileName);
	SF.m_AllowGaps = AllowGaps;

	ProgressStep(0, 1000, "Reading %s", FileName.c_str());
	for (;;)
		{
		unsigned QueryPctDoneX10 = SF.GetPctDoneX10();
		ProgressStep(QueryPctDoneX10, 1000, "Reading %s", FileName.c_str());
		const byte *Seq = SF.GetNextSeq();
		if (Seq == 0)
			break;

		const char *Label = SF.GetLabel();
		unsigned L = SF.GetSeqLength();
		AddSeq(Label, Seq, L);
		}
	ProgressStep(999, 1000, "Reading %s", FileName.c_str());

	SetIsNucleo();

	Progress("%s sequences\n", IntToStr(GetSeqCount()));
	}

void SeqDB::ToFasta(const string &FileName) const
	{
	FILE *f = CreateStdioFile(FileName);
	for (unsigned SeqIndex = 0; SeqIndex < GetSeqCount(); ++SeqIndex)
		ToFasta(f, SeqIndex);
	CloseStdioFile(f);
	}

void SeqDB::SeqToFasta(FILE *f, unsigned SeqIndex, bool WithLabel) const
	{
	if (WithLabel)
		fprintf(f, ">%s\n", GetLabel(SeqIndex));

	const unsigned ROWLEN = 80;

	unsigned L = GetSeqLength(SeqIndex);
	const byte *Seq = GetSeq(SeqIndex);
	unsigned BlockCount = (L + ROWLEN - 1)/ROWLEN;
	for (unsigned BlockIndex = 0; BlockIndex < BlockCount; ++BlockIndex)
		{
		unsigned From = BlockIndex*ROWLEN;
		unsigned To = From + ROWLEN;
		if (To >= L)
			To = L;
		for (unsigned Pos = From; Pos < To; ++Pos)
			fputc(Seq[Pos], f);
		fputc('\n', f);
		}
	}

void SeqDB::ToFasta(FILE *f, unsigned SeqIndex) const
	{
	asserta(SeqIndex < m_SeqCount);
	fprintf(f, ">%s\n", GetLabel(SeqIndex));
	SeqToFasta(f, SeqIndex);
	}

unsigned SeqDB::GetMaxLabelLength() const
	{
	const unsigned SeqCount = GetSeqCount();
	unsigned MaxL = 0;
	for (unsigned Index = 0; Index < SeqCount; ++Index)
		{
		unsigned L = (unsigned) strlen(m_Labels[Index]);
		if (L > MaxL)
			MaxL = L;
		}
	return MaxL;
	}

unsigned SeqDB::GetMaxSeqLength() const
	{
	const unsigned SeqCount = GetSeqCount();
	unsigned MaxL = 0;
	for (unsigned Index = 0; Index < SeqCount; ++Index)
		{
		unsigned L = m_SeqLengths[Index];
		if (L > MaxL)
			MaxL = L;
		}
	return MaxL;
	}

void SeqDB::LogMe() const
	{
	Log("\n");
	const unsigned SeqCount = GetSeqCount();
	Log("SeqDB %u seqs, aligned=%c\n", SeqCount, tof(m_Aligned));
	if (SeqCount == 0)
		return;

	Log("Index             Label  Length  Seq\n");
	Log("-----  ----------------  ------  ---\n");
	for (unsigned Index = 0; Index < SeqCount; ++Index)
		{
		Log("%5u", Index);
		Log("  %16.16s", m_Labels[Index]);
		unsigned L = m_SeqLengths[Index];
		Log("  %6u", L);
		Log("  %*.*s", L, L, m_Seqs[Index]);
		Log("\n");
		}
	}

void SeqDB::GetSeqData(unsigned Id, SeqData &Buffer) const
	{
	asserta(Id < m_SeqCount);
	Buffer.Seq = m_Seqs[Id];
	Buffer.Label = m_Labels[Id];
	Buffer.L = m_SeqLengths[Id];
	Buffer.Index = Id;
	Buffer.ORFParent = 0;
	Buffer.RevComp = false;
	Buffer.Nucleo = IsNucleo();
	}

void SeqDB::SetIsNucleo()
	{
	const unsigned SeqCount = GetSeqCount();
	unsigned N = 0;
	for (unsigned i = 0; i < 100; ++i)
		{
		unsigned SeqIndex = unsigned(rand()%SeqCount);
		const byte *Seq = GetSeq(SeqIndex);
		unsigned L = GetSeqLength(SeqIndex);
		const unsigned Pos = unsigned(rand()%L);
		byte c = Seq[Pos];

		if (g_IsNucleoChar[c])
			++N;
		}
	m_IsNucleo = (N > 80);
	m_IsNucleoSet = true;
	}

unsigned SeqDB::GetTotalLength() const
	{
	const unsigned SeqCount = GetSeqCount();
	unsigned TotalLength = 0;
	for (unsigned Id = 0; Id < SeqCount; ++Id)
		TotalLength += GetSeqLength(Id);
	return TotalLength;
	}

unsigned SeqDB::AddSeq(const char *Label, const byte *Seq, unsigned L)
	{
	StartTimer(AddSeq);
	if (m_SeqCount >= m_Size)
		{
		unsigned NewSize = unsigned(m_Size*1.5) + 1024;
		char **NewLabels = MYALLOC(char *, NewSize, SeqDB);
		byte **NewSeqs = MYALLOC(byte *, NewSize, SeqDB);
		unsigned *NewSeqLengths = MYALLOC(unsigned, NewSize, SeqDB);

		for (unsigned i = 0; i < m_SeqCount; ++i)
			{
			NewLabels[i] = m_Labels[i];
			NewSeqs[i] = m_Seqs[i];
			NewSeqLengths[i] = m_SeqLengths[i];
			}

		MYFREE(m_Labels, m_SeqCount, SeqDB);
		MYFREE(m_Seqs, m_SeqCount, SeqDB);
		MYFREE(m_SeqLengths, m_SeqCount, SeqDB);

		m_Labels = NewLabels;
		m_Seqs = NewSeqs;
		m_SeqLengths = NewSeqLengths;
		m_Size = NewSize;
		}

	unsigned Index = m_SeqCount++;
	m_Seqs[Index] = MYALLOC(byte, L, SeqDB);
	memcpy(m_Seqs[Index], Seq, L);

	unsigned n = strlen(Label) + 1;
	m_Labels[Index] = MYALLOC(char, n, SeqDB);
	memcpy(m_Labels[Index], Label, n);

	if (Index == 0)
		m_Aligned = true;
	else
		m_Aligned = (m_Aligned && L == m_SeqLengths[0]);

	m_SeqLengths[Index] = L;

	EndTimer(AddSeq);
	return Index;
	}

unsigned SeqDB::GetIndex(const char *Label) const
	{
	for (unsigned i = 0; i < m_SeqCount; ++i)
		if (strcmp(Label, m_Labels[i]) == 0)
			return i;
	Die("SeqDB::GetIndex(%s), not found", Label);
	return UINT_MAX;
	}

void SeqDB::MakeLabelToIndex(map<string, unsigned> &LabelToIndex)
	{
	LabelToIndex.clear();
	for (unsigned i = 0; i < m_SeqCount; ++i)
		{
		const string &Label = string(GetLabel(i));
		if (LabelToIndex.find(Label) != LabelToIndex.end())
			Die("Duplicate label: %s", Label.c_str());
		LabelToIndex[Label] = i;
		}
	}
