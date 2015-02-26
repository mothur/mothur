#ifndef mx_h
#define mx_h

#include <list>
#include <limits.h>
#include <math.h>
#include "timing.h"
#include "myutils.h"

const int OPT_LOG = 0x01;
const int OPT_EXP = 0x02;
const int OPT_ZERO_BASED = 0x04;
const float MINUS_INFINITY = -9e9f;
const float UNINIT = -8e8f;

struct SeqData;

template<class T> const char *TypeToStr(T t)
	{
	Die("Unspecialised TypeToStr() called");
	ureturn(0);
	}

template<> inline const char *TypeToStr<unsigned short>(unsigned short f)
	{
	static char s[16];

	sprintf(s, "%12u", f);
	return s;
	}

template<> inline const char *TypeToStr<short>(short f)
	{
	static char s[16];

	sprintf(s, "%12d", f);
	return s;
	}

template<> inline const char *TypeToStr<int>(int f)
	{
	static char s[16];

	sprintf(s, "%5d", f);
	return s;
	}

template<> inline const char *TypeToStr<float>(float f)
	{
	static char s[16];

	if (f == UNINIT)
		sprintf(s, "%12.12s", "?");
	else if (f < MINUS_INFINITY/2)
		sprintf(s, "%12.12s", "*");
	else if (f == 0.0f)
		sprintf(s, "%12.12s", ".");
	else if (f >= -1e5 && f <= 1e5)
		sprintf(s, "%12.5f", f);
	else
		sprintf(s, "%12.4g", f);
	return s;
	}

template<> inline const char *TypeToStr<double>(double f)
	{
	static char s[16];

	if (f < -1e9)
		sprintf(s, "%12.12s", "*");
	else if (f == 0.0f)
		sprintf(s, "%12.12s", ".");
	else if (f >= -1e-5 && f <= 1e5)
		sprintf(s, "%12.5f", f);
	else
		sprintf(s, "%12.4g", f);
	return s;
	}

static inline const char *FloatToStr(float f, string &s)
	{
	s = TypeToStr<float>(f);
	return s.c_str();
	}

template<> inline const char *TypeToStr<char>(char c)
	{
	static char s[2];
	s[0] = c;
	return s;
	}

template<> inline const char *TypeToStr<byte>(byte c)
	{
	static char s[2];
	s[0] = c;
	return s;
	}

template<> inline const char *TypeToStr<bool>(bool tof)
	{
	static char s[2];
	s[0] = tof ? 'T' : 'F';
	return s;
	}

struct SeqDB;

struct MxBase
	{
private:
	MxBase(const MxBase &rhs);
	MxBase &operator=(const MxBase &rhs);

public:
	char m_Name[32];
	char m_Alpha[32];
	char m_Alpha2[32];
	unsigned m_RowCount;
	unsigned m_ColCount;
	unsigned m_AllocatedRowCount;
	unsigned m_AllocatedColCount;
	const SeqDB *m_SeqDB;
	unsigned m_IdA;
	unsigned m_IdB;
	const SeqData *m_SA;
	const SeqData *m_SB;

	static list<MxBase *> *m_Matrices;
	//static MxBase *Get(const string &Name);
	//static float **Getf(const string &Name);
	//static double **Getd(const string &Name);
	//static char **Getc(const string &Name);

	static unsigned m_AllocCount;
	static unsigned m_ZeroAllocCount;
	static unsigned m_GrowAllocCount;
	static double m_TotalBytes;
	static double m_MaxBytes;

	static void OnCtor(MxBase *Mx);
	static void OnDtor(MxBase *Mx);

	MxBase()
		{
		m_AllocatedRowCount = 0;
		m_AllocatedColCount = 0;
		m_RowCount = 0;
		m_ColCount = 0;
		m_IdA = UINT_MAX;
		m_IdB = UINT_MAX;
		m_SeqDB = 0;
		OnCtor(this);
		}
	virtual ~MxBase()
		{
		OnDtor(this);
		}

	virtual unsigned GetTypeSize() const = 0;
	virtual unsigned GetBytes() const = 0;

	void Clear()
		{
		FreeData();
		m_AllocatedRowCount = 0;
		m_AllocatedColCount = 0;
		m_RowCount = 0;
		m_ColCount = 0;
		m_IdA = UINT_MAX;
		m_IdB = UINT_MAX;
		m_SA = 0;
		m_SB = 0;
		}

	bool Empty() const
		{
		return m_RowCount == 0;
		}

	virtual void AllocData(unsigned RowCount, unsigned ColCount) = 0;
	virtual void FreeData() = 0;
	virtual const char *GetAsStr(unsigned i, unsigned j) const = 0;

	void SetAlpha(const char *Alpha)
		{
		unsigned n = sizeof(m_Alpha);
		strncpy(m_Alpha, Alpha, n);
		m_Alpha[n] = 0;
		}

	void Alloc(const char *Name, unsigned RowCount, unsigned ColCount,
	  const SeqDB *DB, unsigned IdA, unsigned IdB,
	  const SeqData *SA, const SeqData *SB);

	void Alloc(const char *Name, unsigned RowCount, unsigned ColCount,
	  const SeqDB *DB = 0, unsigned IdA = UINT_MAX, unsigned IdB = UINT_MAX);

	void Alloc(const char *Name, unsigned RowCount, unsigned ColCount,
	  const SeqData *SA, const SeqData *SB);

	static void LogAll()
		{
		Log("\n");
		if (m_Matrices == 0)
			{
			Log("MxBase::m_Matrices=0\n");
			return;
			}
		Log("\n");
		Log("AllRows  AllCols    Sz        MB  Name\n");
		Log("-------  -------  ----  --------  ----\n");
		double TotalMB = 0;
		for (list<MxBase *>::const_iterator p = m_Matrices->begin();
		  p != m_Matrices->end(); ++p)
			{
			const MxBase *Mx = *p;
			if (Mx == 0)
				continue;
			//if (Mx->m_RowCount != 0 || ShowEmpty)
			//	Mx->LogMe(WithData);
			unsigned ar = Mx->m_AllocatedRowCount;
			if (ar == 0)
				continue;
			unsigned ac = Mx->m_AllocatedColCount;
			unsigned sz = Mx->GetTypeSize();
			double MB = (double) ar*(double) ac*(double) sz/1e6;
			TotalMB += MB;
			Log("%7u  %7u  %4u  %8.2f  %s\n", ar, ac, sz, MB, Mx->m_Name);
			}
		Log("                        --------\n");
		Log("%7.7s  %7.7s  %4.4s  %8.2f\n", "", "", "", TotalMB);
		}

	void LogMe(bool WithData = true, int Opts = 0) const;
	static void LogCounts();
	};

template<class T> struct Mx : public MxBase
	{
// Disable unimplemented stuff
private:
	Mx(Mx &rhs);
	Mx &operator=(Mx &rhs);
	// const Mx &operator=(const Mx &rhs) const;

public:
	T **m_Data;

	Mx()
		{
		m_Data = 0;
		}
	
	~Mx()
		{
		FreeData();
		}

	virtual void AllocData(unsigned RowCount, unsigned ColCount)
		{
		if (opt_logmemgrows)
			Log("MxBase::AllocData(%u,%u) %s bytes, Name=%s\n",
			  RowCount, ColCount, IntToStr(GetBytes()), m_Name);
		// m_Data = myalloc<T *>(RowCount);
		m_Data = MYALLOC(T *, RowCount, Mx);
		for (unsigned i = 0; i < RowCount; ++i)
			// m_Data[i] = myalloc<T>(ColCount);
			m_Data[i] = MYALLOC(T, ColCount, Mx);
		AddBytes("Mx_AllocData", RowCount*sizeof(T *) + RowCount*ColCount*sizeof(T));

		m_AllocatedRowCount = RowCount;
		m_AllocatedColCount = ColCount;
		}

	virtual void FreeData()
		{
		for (unsigned i = 0; i < m_AllocatedRowCount; ++i)
			MYFREE(m_Data[i], m_AllocatedColCount, Mx);
		MYFREE(m_Data, m_AllocatedRowCount, Mx);
		SubBytes("Mx_AllocData",
		  m_AllocatedRowCount*sizeof(T *) + m_AllocatedRowCount*m_AllocatedColCount*sizeof(T));

		m_Data = 0;
		m_RowCount = 0;
		m_ColCount = 0;
		m_AllocatedRowCount = 0;
		m_AllocatedColCount = 0;
		}

	T **GetData()
		{
		return (T **) m_Data;
		}

	T Get(unsigned i, unsigned j) const
		{
		assert(i < m_RowCount);
		assert(j < m_ColCount);
		return m_Data[i][j];
		}

	void Put(unsigned i, unsigned j, T x) const
		{
		assert(i < m_RowCount);
		assert(j < m_ColCount);
		m_Data[i][j] = x;
		}

	T GetOffDiagAvgs(vector<T> &Avgs) const
		{
		if (m_RowCount != m_ColCount)
			Die("GetOffDiagAvgs, not symmetrical");
		Avgs.clear();
		T Total = T(0);
		for (unsigned i = 0; i < m_RowCount; ++i)
			{
			T Sum = T(0);
			for (unsigned j = 0; j < m_ColCount; ++j)
				{
				if (j == i)
					continue;
				Sum += m_Data[i][j];
				}
			T Avg = Sum/(m_RowCount-1);
			Total += Avg;
			Avgs.push_back(Avg);
			}
		return m_RowCount == 0 ? T(0) : Total/m_RowCount;
		}

	unsigned GetTypeSize() const
		{
		return sizeof(T);
		}

	virtual unsigned GetBytes() const
		{
		return m_AllocatedRowCount*m_AllocatedColCount*GetTypeSize() +
		  m_AllocatedRowCount*sizeof(T *);
		}

	const char *GetAsStr(unsigned i, unsigned j) const
		{
		return TypeToStr<T>(Get(i, j));
		}

	const T *const *const GetData() const
		{
		return (const T *const *) m_Data;
		}

	void Copy(const Mx<T> &rhs)
		{
		Alloc("Copy", rhs.m_RowCount, rhs.m_ColCount, rhs.m_SeqDB, rhs.m_IdA, rhs.m_IdB);
		const T * const *Data = rhs.GetData();
		for (unsigned i = 0; i < m_RowCount; ++i)
			for (unsigned j = 0; j < m_ColCount; ++j)
				m_Data[i][j] = Data[i][j];
		}

	void Assign(T v)
		{
		for (unsigned i = 0; i < m_RowCount; ++i)
			for (unsigned j = 0; j < m_ColCount; ++j)
				m_Data[i][j] = v;
		}

	bool Eq(const Mx &rhs, bool Bwd = false) const
		{
		if (rhs.m_ColCount != m_ColCount)
			return false;
		if (rhs.m_RowCount != m_RowCount)
			return false;
		const T * const*d = rhs.GetData();
		int i1 = Bwd ? m_RowCount : 0;
		int j1 = Bwd ? m_ColCount : 0;
		int i2 = Bwd ? -1 : m_RowCount;
		int j2 = Bwd ? -1 : m_ColCount;
		for (int i = i1; i != i2; Bwd ? --i : ++i)
			for (int j = j1; j != j2; Bwd ? --j : ++j)
				{
				float x = m_Data[i][j];
				float y = d[i][j];
				if (x < -1e10 && y < -1e10)
					continue;
				if (!feq(x, y))
					{
					Warning("%s[%d][%d] = %g, %s = %g",
					  m_Name, i, j, x, rhs.m_Name, y);
					return false;
					}
				}
		return true;
		}

	bool EqMask(const Mx &rhs, const Mx<bool> &Mask) const
		{
		if (rhs.m_ColCount != m_ColCount)
			return false;
		if (rhs.m_RowCount != m_RowCount)
			return false;

		if (Mask.m_ColCount != m_ColCount)
			return false;
		if (Mask.m_RowCount != m_RowCount)
			return false;

		const T * const*d = rhs.GetData();
		bool Bwd = false;
		int i1 = Bwd ? m_RowCount : 0;
		int j1 = Bwd ? m_ColCount : 0;
		int i2 = Bwd ? -1 : m_RowCount;
		int j2 = Bwd ? -1 : m_ColCount;
		for (int i = i1; i != i2; Bwd ? --i : ++i)
			for (int j = j1; j != j2; Bwd ? --j : ++j)
				{
				if (!Mask.m_Data[i][j])
					continue;
				float x = m_Data[i][j];
				float y = d[i][j];
				if (x < -1e10 && y < -1e10)
					continue;
				if (!feq(x, y))
					{
					Warning("%s[%d][%d] = %g, %s = %g",
					  m_Name, i, j, x, rhs.m_Name, y);
					return false;
					}
				}
		return true;
		}

	void Init(T v)
		{
		for (unsigned i = 0; i < m_RowCount; ++i)
			for (unsigned j = 0; j < m_ColCount; ++j)
				m_Data[i][j] = v;
		}
	};

void WriteMx(const string &Name, Mx<float> &Mxf);

template<class T> void ReserveMx(Mx<T> &Mxf, unsigned N = UINT_MAX)
	{
	if (Mxf.m_AllocatedRowCount > 0)
		return;
	extern unsigned g_MaxInputSeqLength;
	if (N == UINT_MAX)
		N = g_MaxInputSeqLength+1;
	Mxf.Alloc("(Reserved)", N, N);
	}

#endif // mx_h
