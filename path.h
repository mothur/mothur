#ifndef path_h
#define path_h

struct PathData
	{
private:
	PathData(PathData &);
	PathData &operator=(PathData &);

public:
	char *Start;
	char *Front;
	char *Back;
	unsigned Bytes;

public:
	PathData()
		{
		Clear(true);
		}
	~PathData()
		{
		Free();
		}
	void Free();
	void Alloc(unsigned MaxLen);
	void Clear(bool ctor = false)
		{
		Start = 0;
		if (ctor)
			{
			Front = 0;
			Back = 0;
			Bytes = 0;
			}
		else
			Free();
		}
	void Copy(const PathData &rhs);
	void FromStr(const char *PathStr);
	void Reverse()
		{
		asserta(Start != 0);
		unsigned L = (unsigned) strlen(Start);
		for (unsigned k = 0; k < L/2; ++k)
			{
			char c = Start[k];
			Start[k] = Start[L-k-1];
			Start[L-k-1] = c;
			}
		}
	void SetEmpty()
		{
		Start = 0;
		}

	bool IsEmpty() const
		{
		return Start == 0;
		}
	};

#endif // path_h
