#ifndef sfasta_h
#define sfasta_h

#include "myutils.h"
#include "seq.h"

typedef void (*ON_START_XSEQ)(const SeqData &SD);
typedef void (*ON_END_XSEQ)(const SeqData &SD);

// Sequential reader for FASTA file format.
// Serves sequences in file order to save memory.
// Caches biggish chunks to compromise memory vs. speed.
class SFasta
	{
public:
	string m_FileName;
	FILE *m_File;
	bool m_AllowGaps;

	off_t m_FileSize;

// Position to start next read
	off_t m_FilePos;

// Cached data.
	byte *m_Buffer;

// Bytes allocated to m_Buffer
	unsigned m_BufferSize;

// Current position in buffer, normally points to '>'
	unsigned m_BufferOffset;

// File data in buffer <= m_BufferSize
	unsigned m_BufferBytes;

// Current label
// Points into m_Buffer, not a separate buffer.
	char *m_Label;

// Current sequence length
	unsigned m_SeqLength;

// Current seq index
	unsigned m_SeqIndex;

	unsigned m_ShortestLength;
	unsigned m_LongestLength;
	unsigned m_TooShortCount;
	unsigned m_TooLongCount;
	unsigned m_TooPolyCount;

private:
	bool m_IsNucleoSet;
	bool m_IsNucleo;

public:
	SFasta();
	~SFasta();

	void Clear();
	void Open(const string &FileName);
	void Rewind();
	bool SetIsNucleo();
	bool GetIsNucleo() const { asserta(m_IsNucleoSet); return m_IsNucleo; };

// Get next sequence.
// Returns zero on end-of-file
	const byte *GetNextSeq();

// Get next sequence as SeqData object, return false on end-of-file.
	bool GetNextSD(SeqData &SD);

// Length of most recent sequence returned by GetNextSeq().
	unsigned GetSeqLength() const { return m_SeqLength; }

// Label of most recent sequence returned by GetNextSeq().
	const char *GetLabel() const { return m_Label; }

// Index of most recent sequence returned by GetNextSeq().
	unsigned GetSeqIndex() const { return m_SeqIndex; }

	unsigned GetPctDoneX10() const;
	double GetPctDone() const;

	void LogMe() const;

private:
	void FillCache();
	const byte *GetNextSeqLo();
	};

#endif // sfasta_h
