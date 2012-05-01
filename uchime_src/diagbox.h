#ifndef diagbox_h
#define diagbox_h

struct DiagBox;

void GetDiagBox(unsigned LA, unsigned LB, unsigned DiagLo, unsigned DiagHi, DiagBox &Box);
void GetDiagRange(unsigned LA, unsigned LB, unsigned d,
  unsigned &mini, unsigned &minj, unsigned &maxi, unsigned &maxj);
void GetDiagLoHi(unsigned LA, unsigned LB, const char *Path,
  unsigned &dlo, unsigned &dhi);

struct DiagBox
	{
	DiagBox()
		{
		}

	DiagBox(unsigned LA_, unsigned LB_, unsigned DiagLo, unsigned DiagHi)
		{
		//GetDiagBox(LA, LB, DiagLo, DiagHi, *this);
		//Validate();
		Init(LA_, LB_, DiagLo, DiagHi);
		}

	void Init(unsigned LA_, unsigned LB_, unsigned DiagLo, unsigned DiagHi)
		{
		GetDiagBox(LA_, LB_, DiagLo, DiagHi, *this);
		Validate();
		}

	unsigned LA;
	unsigned LB;

	unsigned dlo;
	unsigned dhi;

	unsigned dlo_mini;
	unsigned dlo_minj;

	unsigned dlo_maxi;
	unsigned dlo_maxj;

	unsigned dhi_mini;
	unsigned dhi_minj;

	unsigned dhi_maxi;
	unsigned dhi_maxj;

	unsigned GetDiag(unsigned i, unsigned j) const
		{
		return LA - i + j;
		}

// i, j are positions 0..LA-1, 0..LB-1.
	bool InBox(unsigned i, unsigned j) const
		{
		unsigned d = GetDiag(i, j);
		return d >= dlo && d <= dhi;
		}

/***
i, j are 0-based prefix lengths 0..LA, 0..LB.

A full path is in the box iff all match pairs are in the box.

A partial path that aligns a prefix of A to a prefix of B as
in D.P.) is in the box iff it is is the prefix of at least
one full path that is in the box.

A D.P. matrix entry X[i][j] is in the box iff there is at
least one full path aligning the first i letters of A and
the first j letters of B ending in a column of type X, i.e.
if there exists a partial path in the box that ends in X.

Assume terminals appear in all paths, and DI/ID forbidden.

Intuitively seems that by these definitions D is in box iff
DM or MD is in box, I is in box iff IM or MI is in box.
Don't have proof..
***/
	bool InBoxDPM(unsigned i, unsigned j) const
		{
	// Special case for M[0][0]
		if (i == 0 && j == 0)
			return true;
		if (i == 0 || j == 0)
			return false;
		unsigned d = GetDiag(i-1, j-1);
		return d >= dlo && d <= dhi;
		}

	bool InBoxDPD(unsigned i, unsigned j) const
		{
		bool MD = i == 0 ? false : InBoxDPM(i-1, j);
		bool DM = (i == LA || j == LB) ? false : InBoxDPM(i+1, j+1);
		return MD || DM;
		}

	bool InBoxDPI(unsigned i, unsigned j) const
		{
		bool MI = j == 0 ? false : InBoxDPM(i, j-1);
		bool IM = (i == LA || j == LB) ? false : InBoxDPM(i+1, j+1);
		return MI || IM;
		}

	// d = LA - i + j = 1 .. LA+LB-1
	void Validate() const
		{
		asserta(dlo <= dhi);
		asserta(dlo >= GetDiag(LA-1, 0));
		asserta(dhi <= GetDiag(0, LB-1));

		asserta(GetDiag(dlo_mini, dlo_minj) == dlo);
		asserta(GetDiag(dlo_maxi, dlo_maxj) == dlo);
		asserta(GetDiag(dhi_mini, dhi_minj) == dhi);
		asserta(GetDiag(dhi_maxi, dhi_maxj) == dhi);

		asserta(dlo_mini >= dhi_mini);
		asserta(dlo_minj <= dhi_minj);
		asserta(dlo_maxi >= dhi_maxi);
		asserta(dlo_maxj <= dhi_maxj);
		}

	unsigned GetMini() const
		{
		return dhi_mini;
		}

	unsigned GetMaxi() const
		{
		return dlo_maxi;
		}

	unsigned GetMinj() const
		{
		return dlo_minj;
		}

	unsigned GetMaxj() const
		{
		return dhi_maxj;
		}
/***
	i = 0..LA-1
	j = 0..LB-1
	d = LA - i + j = 1 .. LA+LB-1
	j = d - LA + i
	i = LA - d + j
***/
	void GetRange_j(unsigned i, unsigned &Startj, unsigned &Endj) const
		{
	// j = d - LA + i
		if (dlo + i >= LA)
			Startj = dlo + i - LA;
		else
			Startj = 0;

		if (Startj >= LB)
			Startj = LB - 1;

		if (dhi + i + 1 >= LA)
			Endj = dhi + i + 1 - LA;
		else
			Endj = 0;

		if (Endj > LB)
			Endj = LB;

		asserta(Endj >= Startj);
		}

	void LogMe() const
		{
		Log("LA=%u LB=%d dlo(%u): (%u,%u)-(%u,%u) dhi(%u): (%u,%u)-(%u,%u) i=[%u-%u] j=[%u-%u]\n",
		  LA, LB,
		  dlo,
		  dlo_mini, dlo_minj,
		  dlo_maxi, dlo_maxj,
		  dhi,
		  dhi_mini, dhi_minj,
		  dhi_maxi, dhi_maxj,
		  GetMini(), GetMaxi(),
		  GetMinj(), GetMaxj());
		}
	};

typedef const char *(*NWDIAG)(const byte *A, unsigned LA, const byte *B, unsigned LB,
  unsigned DiagLo, unsigned DiagHi, bool LeftTerm, bool RightTerm);

const char *NWBandWrap(NWDIAG NW, const byte *A, unsigned LA, const byte *B, unsigned LB,
  unsigned DiagLo, unsigned DiagHi, bool LeftTerm, bool RightTerm);

#endif // diagbox_h
