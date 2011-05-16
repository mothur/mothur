#ifndef seq_h
#define seq_h

struct ORFData;

struct SeqData
	{
	const char *Label;
	const byte *Seq;
	unsigned L;
	unsigned Index;

// RevComp means that SeqData.Seq is reverse-complemented relative
// to the sequence in the input file (query or db). Coordinates in
// a hit (e.g., AlnData) will be relative to SeqData.Seq, so both
// the sequence and the coordinates should be r.c.'d for output.
	bool RevComp;
	bool Nucleo;
	const ORFData *ORFParent;

	SeqData()
		{
		Clear();
		}

	void Clear()
		{
		Label = 0;
		Seq = 0;
		L = 0;
		Index = UINT_MAX;
		RevComp = false;
		Nucleo = false;
		ORFParent = 0;
		}
	};

#endif // seq_h
