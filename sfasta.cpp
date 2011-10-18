#include "sfasta.h"
#include "orf.h"
#include "alpha.h"
#include "timing.h"

static inline bool isgap(byte c)
	{
	return c == '-' || c == '.';
	}

const unsigned BufferSize = 16*1024*1024;

static unsigned GetMaxPoly(const byte *Seq, unsigned L)
	{
	byte CurrChar = Seq[0];
	unsigned Start = 0;
	unsigned MaxLen = 1;
	for (unsigned i = 1; i < L; ++i)
		{
		char c = Seq[i];
		if (c != CurrChar || i+1 == L)
			{
			unsigned Len = i - Start;
			if (Len > MaxLen)
				MaxLen = Len;
			CurrChar = c;
			Start = i;
			}
		}
	return MaxLen;
	}

SFasta::SFasta()
	{
	m_FileName = "";
	m_File = 0;
	m_Buffer = 0;
	m_BufferSize = 0;
	m_BufferOffset = 0;
	m_BufferBytes = 0;
	m_FilePos = 0;
	m_FileSize = 0;
	m_Label = 0;
	m_SeqLength = 0;
	m_TooShortCount = 0;
	m_TooLongCount = 0;
	m_ShortestLength = 0;
	m_LongestLength = 0;
	m_IsNucleo = false;
	m_IsNucleoSet = false;
	}

SFasta::~SFasta()
	{
	Clear();
	}

void SFasta::Clear()
	{
	MYFREE(m_Buffer, m_BufferSize, SFasta);
	if (m_File != 0)
		CloseStdioFile(m_File);

	m_FileName = "";
	m_File = 0;
	m_Buffer = 0;
	m_BufferSize = 0;
	m_BufferOffset = 0;
	m_BufferBytes = 0;
	m_FilePos = 0;
	m_FileSize = 0;
	m_Label = 0;
	m_SeqLength = 0;
	m_SeqIndex = UINT_MAX;
	m_AllowGaps = false;
	m_IsNucleo = false;
	m_IsNucleoSet = false;
	m_TooShortCount = 0;
	m_TooLongCount = 0;
	m_ShortestLength = 0;
	m_LongestLength = 0;
	m_TooPolyCount = 0;
	}

void SFasta::LogMe() const
	{
	Log("\n");
	Log("SFasta::LogMe()\n");
	Log("FileName=%s\n", m_FileName.c_str());
	Log("FileSize=%u\n", (unsigned) m_FileSize);
	Log("FilePos=%u\n", (unsigned) m_FilePos);
	Log("BufferSize=%u\n", m_BufferSize);
	Log("BufferPos=%u\n", m_BufferOffset);
	Log("BufferBytes=%u\n", m_BufferBytes);
	if (m_Label == 0)
		Log("Label=NULL\n");
	else
		Log("Label=%s\n", m_Label);
	Log("SeqLength=%u\n", m_SeqLength);
	}

const byte *SFasta::GetNextSeq()
	{
	for (;;)
		{
		const byte *Seq = GetNextSeqLo();
		if (Seq == 0)
			{
			if (m_TooShortCount > 0)
				Warning("%u short sequences (--minlen %u, shortest %u) discarded from %s",
				  m_TooShortCount, opt_minlen, m_ShortestLength, m_FileName.c_str());
			if (m_TooLongCount > 0)
				Warning("%u long sequences (--maxlen %u, longest %u) discarded from %s",
				  m_TooLongCount, opt_maxlen, m_LongestLength, m_FileName.c_str());
			if (m_TooPolyCount > 0)
				Warning("%u sequences with long homopolymers discarded (--maxpoly %u)",
				  m_TooPolyCount, opt_maxpoly);
			return 0;
			}
		if (m_SeqLength < opt_minlen)
			{
			++m_TooShortCount;
			if (m_ShortestLength == 0 || m_SeqLength < m_ShortestLength)
				m_ShortestLength = m_SeqLength;
			continue;
			}
		if (m_SeqLength > opt_maxlen && opt_maxlen != 0)
			{
			if (m_LongestLength == 0 || m_SeqLength > m_LongestLength)
				m_LongestLength = m_SeqLength;
			++m_TooLongCount;
			continue;
			}
		return Seq;
		}
	}

const byte *SFasta::GetNextSeqLo()
	{
// End of cache?
	if (m_BufferOffset == m_BufferBytes)
		{
	// End of file?
		if (m_FilePos == m_FileSize)
			return 0;
		FillCache();
		}

	StartTimer(SF_GetNextSeq);
	asserta(m_Buffer[m_BufferOffset] == '>');
	m_Label = (char *) (m_Buffer + m_BufferOffset + 1);
	
//// Scan to end-of-line.
//// Use dubious library function strchr() in the hope
//// that it uses fast machine code.
//	byte *ptr = (byte *) strchr(m_Label, '\n');
//	asserta(ptr != 0);
//	*ptr = 0;

	byte *ptr = 0;
	for (unsigned i = m_BufferOffset; i < m_BufferSize; ++i)
		{
		char c = m_Buffer[i];
		if (c == '\n' || c == '\r')
			{
			ptr = m_Buffer + i;
			break;
			}
		}
	asserta(ptr != 0);

	if (opt_trunclabels)
		{
		for (char *p = m_Label; *p; ++p)
			if (isspace(*p))
				{
				*p = 0;
				break;
				}
		}
	else
		{
		for (char *p = m_Label; *p; ++p)
			{
			if (*p == '\t')
				*p = ' ';
			else if (*p == '\r' || *p == '\n')
				{
				*p = 0;
				char NextChar = *(p+1);
				if (NextChar == '\r' || NextChar == '\n')
					++p;
				break;
				}
			}
		}

// ptr points to end-of-line.
// Move to start of sequence data.
	byte *Seq = ++ptr;

// Delete white space in-place
	byte *To = ptr;
	m_BufferOffset = (unsigned) (ptr - m_Buffer);
	while (m_BufferOffset < m_BufferBytes)
		{
		byte c = m_Buffer[m_BufferOffset];
		if (c == '>')
			{
			char prevc = '\n';
			if (m_BufferOffset > 0)
				prevc = m_Buffer[m_BufferOffset-1];
			if (prevc == '\n' || prevc == '\r')
				break;
			}
		++m_BufferOffset;
		if (isalpha(c) || (isgap(c) && m_AllowGaps))
			*To++ = c;
		else if (c == '\n' || c == '\r')
			continue;
		else
			{
			const char *Label = (m_Label == 0 ? "" : m_Label);
			static bool WarningDone = false;
			if (!WarningDone)
				{
				if (isgap(c))
					Warning("Ignoring gaps in FASTA file '%s'",
					  m_FileName.c_str());
				else if (isprint(c))
					Warning("Invalid FASTA file '%s', non-letter '%c' in sequence >%s",
					  m_FileName.c_str(), c, Label);
				else
					Warning("Invalid FASTA file '%s', non-printing byte (hex %02x) in sequence >%s",
					  m_FileName.c_str(), c, Label);
				WarningDone = true;
				}
			continue;
			}
		}
	m_SeqLength = unsigned(To - Seq);

	if (m_SeqIndex == UINT_MAX)
		m_SeqIndex = 0;
	else
		++m_SeqIndex;

	EndTimer(SF_GetNextSeq);
	return Seq;
	}

void SFasta::Open(const string &FileName)
	{
	Clear();
	m_FileName = FileName;
	m_File = OpenStdioFile(FileName);
	m_BufferSize = BufferSize;
	//m_Buffer = myalloc<byte>(m_BufferSize);
	m_Buffer = MYALLOC(byte, m_BufferSize, SFasta);
	m_FileSize = GetStdioFileSize(m_File);
	}

void SFasta::Rewind()
	{
	m_BufferOffset = 0;
	m_BufferBytes = 0;
	m_FilePos = 0;
	}

bool SFasta::SetIsNucleo()
	{
	if (m_FilePos != 0)
		Die("SFasta::IsNucleo, not at BOF");

	unsigned LetterCount = 0;
	unsigned NucleoLetterCount = 0;
	for (;;)
		{
		const byte *Seq = GetNextSeq();
		if (Seq == 0)
			break;
		unsigned L = GetSeqLength();
		for (unsigned i = 0; i < L; ++i)
			if (g_IsNucleoChar[Seq[i]])
				++NucleoLetterCount;
		LetterCount += L;
		if (LetterCount > 256)
			break;
		}
	Rewind();
	if (LetterCount == 0)
		{
		m_IsNucleoSet = true;
		m_IsNucleo = true;
		return true;
		}

// Nucleo if more than 90% nucleo letters AGCTUN
	m_IsNucleo = double(NucleoLetterCount)/LetterCount > 0.9;
	m_IsNucleoSet = true;
	return m_IsNucleo;
	}

void SFasta::FillCache()
	{
	StartTimer(SF_FillCache);
	asserta(m_FilePos < m_FileSize);

// off_t may be larger type than unsigned, e.g. 64- vs. 32-bit.
	off_t otBytesToRead = m_FileSize - m_FilePos;

	bool FinalBuffer = true;
	if (otBytesToRead > (off_t) m_BufferSize)
		{
		FinalBuffer = false;
		otBytesToRead = m_BufferSize;
		}

	unsigned BytesToRead = unsigned(otBytesToRead);
	asserta(BytesToRead > 0);
	asserta(BytesToRead <= m_BufferSize);

	SetStdioFilePos(m_File, m_FilePos);
	ReadStdioFile(m_File, m_Buffer, BytesToRead);
	if (m_Buffer[0] != '>')
		{
		if (m_FilePos == 0)
			Die("Input is not FASTA file");
		else
			Die("SFasta::FillCache() failed, expected '>'");
		}

	m_BufferOffset = 0;

// If last buffer in file, done
	if (FinalBuffer)
		{
		m_BufferBytes = BytesToRead;
		m_FilePos += BytesToRead;
		EndTimer(SF_FillCache);
		return;
		}

// If not last buffer, truncate any partial sequence
// at end of buffer. Search backwards to find last '>'.
	byte *ptr = m_Buffer + BytesToRead - 1;
	while (ptr > m_Buffer)
		{
		if (ptr[0] == '>' && (ptr[-1] == '\n' || ptr[-1] == '\r'))
			break;
		--ptr;
		}

	if (ptr == m_Buffer)
		{
		LogMe();
		if (*ptr != '>')
			{
	// No '>' found.
	// This might techincally be legal FASTA if the entire
	// buffer is white space, but strange if not the last buffer
	// in the file, so quit anyway.
			Die("Failed to find '>' (pos=%u, bytes=%u)",
			  (unsigned) m_FilePos, BytesToRead);
			}
		else
			{
	// Entire buffer is one sequence which may be truncated.
			Die("Sequence too long (pos=%u, bytes=%u)",
			  (unsigned) m_FilePos, BytesToRead);
			}
		}

	asserta(*ptr == '>');

	m_BufferBytes = unsigned(ptr - m_Buffer);
	m_FilePos += m_BufferBytes;

	EndTimer(SF_FillCache);
	}

unsigned SFasta::GetPctDoneX10() const
	{
	if (m_FilePos == 0 || m_FileSize == 0)
		return 0;

	assert(m_FilePos >= (off_t) m_BufferBytes);
	off_t BufferStart = m_FilePos - m_BufferBytes;
	off_t BufferPos = BufferStart + m_BufferOffset;

	unsigned iPctX10 = unsigned(10.0*double(BufferPos)*100.0/double(m_FileSize));
	if (iPctX10 == 0)
		return 1;
	if (iPctX10 >= 999)
		return 998;
	return iPctX10;
	}

double SFasta::GetPctDone() const
	{
	if (m_FilePos == 0 || m_FileSize == 0)
		return 0;

	assert(m_FilePos >= (off_t) m_BufferBytes);
	off_t BufferStart = m_FilePos - m_BufferBytes;
	off_t BufferPos = BufferStart + m_BufferOffset;

	return double(BufferPos)*100.0/double(m_FileSize);
	}

bool SFasta::GetNextSD(SeqData &SD)
	{
	SD.Seq = GetNextSeq();
	if (SD.Seq == 0)
		return false;

	SD.Label = GetLabel();
	SD.L = GetSeqLength();
	SD.Index = GetSeqIndex();
	SD.ORFParent = 0;
	SD.Nucleo = GetIsNucleo();
	SD.RevComp = false;

	return true;
	}

#if	TEST
void TestSFasta()
	{
	SFasta SF;
	SF.Open(opt_input);

	if (opt_verbose)
		{
		Log("  Index   Length  Label\n");
		Log("-------  -------  -----\n");
		}

	unsigned Index = 0;
	unsigned SeqCount = 0;
	double LetterCount = 0.0;
	ProgressStep(0, 1000, "Reading");
	for (;;)
		{
		const byte *Seq = SF.GetNextSeq();
		if (Seq == 0)
			break;
		ProgressStep(SF.GetPctDoneX10(), 1000, "Reading");
		const char *Label = SF.GetLabel();
		unsigned L = SF.GetSeqLength();
		++SeqCount;
		LetterCount += L;

		if (opt_verbose)
			{
			Log(">%7u  %7u  '%s'\n", Index, L, Label);
			Log("+%7.7s  %7.7s  \"%*.*s\"\n", "", "", L, L, Seq);
			}

		++Index;
		}
	ProgressStep(999, 1000, "Reading");

	Progress("%u seqs, %s letters\n", SeqCount, FloatToStr(LetterCount));
	Log("%u seqs, %s letters\n", SeqCount, FloatToStr(LetterCount));
	}
#endif // TEST
