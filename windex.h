//uchime by Robert C. Edgar http://drive5.com/uchime This code is donated to the public domain.

#ifndef windex_h
#define windex_h

class SFasta;
struct SeqDB;

typedef uint32 word_t;
typedef uint16 wordcount_t;
typedef uint32 arrsize_t;
typedef uint16 seqcountperword_t;
typedef uint32 seqindex_t;
typedef uint16 commonwordcount_t;

const uint32 WindexFileHdr_Magic1 = 0x312DE41;
const uint32 WindexFileHdr_Magic2 = 0x312DE42;
const uint32 WindexFileHdr_Magic3 = 0x312DE43;
const uint32 WindexFileHdr_Magic4 = 0x312DE44;

struct WindexFileHdr
	{
	uint32 Magic1;
	uint32 IsNucleo;
	uint32 WordLength;
	uint32 Magic2;
	};

class Windex
	{
public:
	bool m_Nucleo;
	bool m_RedAlpha;
	unsigned m_WordLength;
	unsigned m_AlphaSize;
	unsigned m_WordCount;
	unsigned m_Hi;
	unsigned m_CapacityInc;
	arrsize_t *m_Capacities;
	arrsize_t *m_Sizes;
	float *m_WordScores;
	seqindex_t **m_SeedIndexes;
	byte *m_UniqueCounts;
	unsigned m_CharToLetter[256];

public:
	Windex();
	void ToFile(const string &FileName) const;
	void FromFile(const string &FileName);
	void FromSFasta(SFasta &SF);
	void FromSeqDB(const SeqDB &DB);
	void Clear(bool ctor = false);
	void AddWords(unsigned SeqIndex, const word_t *Words, unsigned N);
	void Init(bool Nucleo, unsigned WordLength);
	void Init2(bool Nucleo, unsigned TableSize);
	void InitRed(unsigned WordLength);
	void InitWordScores(const float *const *SubstMx);
	void Reset();
	void LogMe() const;
	unsigned LogMemSize() const;
	void LogWordStats(unsigned TopWords = 10) const;
	const char *WordToStr(word_t Word) const;
	word_t SeqToWord(const byte *Seq) const;
	unsigned SeqToWords(const byte *Seq, unsigned L, word_t *Words) const;
	unsigned SeqToWordsStep(unsigned Step, const byte *Seq, unsigned L, word_t *Words) const;
	unsigned WordsToCounts(const word_t *Words, unsigned N,
	  word_t *UniqueWords, seqcountperword_t *Counts) const;
	unsigned GetUniqueWords(const word_t *Words, unsigned N,
	  word_t *UniqueWords) const;
	void LogSizeHisto() const;
	};

#endif // windex_h
