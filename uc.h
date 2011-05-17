//uchime by Robert C. Edgar http://drive5.com/uchime This code is donated to the public domain.

#ifndef uc_h
#define uc_h

#include "seqdb.h"
#include "seq.h"
#include "path.h"

struct AlnData;

int uchime_main(int, char**);  

class UCFile
	{
public:
	FILE *m_File;
	byte *m_Data;
	vector<char> m_RecTypes;
	vector<float> m_PctIds;
	vector<const char *> m_Labels;
	vector<const char *> m_SeedLabels;
	vector<unsigned> m_SeedIndexes;
	vector<const char *> m_CompressedPaths;
	vector<unsigned> m_SeqLengths;
	vector<unsigned> m_SortOrder;
	vector<char> m_Strands;
	vector<unsigned> m_Los;
	vector<unsigned> m_SeedLos;

public:
	/* some function prototypes */
	
		
	UCFile();
	void Clear(bool ctor = false);
	void Close();
	void FromFile(const string &FileName);
	void FromClstr(const string &FileName);
	void ToFile(const string &FileName);
	unsigned GetRecordCount() const;
	void LogMe() const;
	void ToClstr(const string &FileName);
	void ToFasta(const string &FileName, const SeqDB &Input, bool Reformat);
	void Create(const string &FileName);
	void Sort();
	void Flush() const;

	void WriteNotMatched(unsigned L, const char *Label) const;
	void WriteLibSeed(unsigned SeedIndex, unsigned L, const char *Label) const;
	void WriteNewSeed(unsigned SeedIndex, unsigned L, const char *Label) const;
	void WriteHit(const SeqData &SA, const SeqData &SB, double FractId,
	  const PathData &PD) const;
	void WriteReject(const SeqData &SA, const SeqData &SB, double FractId,
	  const char *Path) const;
	void WriteHit(unsigned SeedIndex, unsigned L, double PctId,
	  const char *CompressedPath, char Strand, unsigned Lo, unsigned SeedLo,
	  const char *Label, const char *SeedLabel) const;
	void WriteHit(const AlnData &AD);
	void WriteLibCluster(unsigned SeedIndex, unsigned Size, double AvgId,
	  const char *Label) const;
	void WriteNewCluster(unsigned SeedIndex, unsigned Size, double AvgId,
	  const char *Label) const;
	void WriteSeqX(FILE *f, const byte *Seq, unsigned L, const char *CompressedPath) const;
	};

#endif // uc_h
