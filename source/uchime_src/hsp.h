#ifndef hsp_h
#define hsp_h	1

struct HSPData
	{
	unsigned Loi;
	unsigned Loj;
	unsigned Leni;
	unsigned Lenj;
	float Score;
	unsigned User;

	unsigned GetLength() const
		{
		if (Leni != Lenj)
			Die("HSP::GetLength(): Leni %u, Lenj %u, Loi %u, Loj %u, Score %.1f",
			  Leni, Lenj, Loi, Loj, Score);

		return Leni;
		}

	unsigned GetHii() const
		{
		assert(Leni > 0);
		return Loi + Leni - 1;
		}

	unsigned GetHij() const
		{
		assert(Lenj > 0);
		return Loj + Lenj - 1;
		}

	bool LeftA() const
		{
		return Loi == 0;
		}

	bool LeftB() const
		{
		return Loj == 0;
		}

	bool RightA(unsigned LA) const
		{
		return Loi + Leni == LA;
		}

	bool RightB(unsigned LB) const
		{
		return Loj + Lenj == LB;
		}

	unsigned GetIdCount(const byte *A, const byte *B) const
		{
		unsigned Count = 0;
		unsigned K = GetLength();
		for (unsigned k = 0; k < K; ++k)
			{
			byte a = A[Loi+k];
			byte b = B[Loj+k];
			if (toupper(a) == toupper(b))
				Count++;
			}
		return Count;
		}

	double OverlapFract(const HSPData &HSP) const
		{
		if (Leni == 0 || Lenj == 0)
			return 0.0;

		unsigned MaxLoi = max(Loi, HSP.Loi);
		unsigned MaxLoj = max(Loj, HSP.Loj);
		unsigned MinHii = min(GetHii(), HSP.GetHii());
		unsigned MinHij = min(GetHij(), HSP.GetHij());

		unsigned Ovi = (MinHii < MaxLoi) ? 0 : MinHii - MaxLoi;
		unsigned Ovj = (MinHij < MaxLoj) ? 0 : MinHij - MaxLoj;

		asserta(Ovi <= Leni && Ovj <= Lenj);
		return double(Ovi*Ovj)/double(Leni*Lenj);
		}

	bool operator<(const HSPData &rhs) const
		{
		return Loi < rhs.Loi;
		}

	void LogMe() const
		{
		Log("Loi=%u Loj=%u Li=%u Lj=%u Score=%.1f\n", Loi, Loj, Leni, Lenj, Score);
		}

	void LogMe2() const
		{
		Log("(%u-%u,%u-%u/%.1f)", Loi, GetHii(), Loj, GetHij(), Score);
		}
	};

// Bendpoint
struct BPData
	{
	unsigned Pos;
	bool IsLo;
	unsigned Index;

	void LogMe() const
		{
		Log("BP%s Pos %u Ix %u", (IsLo ? "lo" : "hi"), Pos, Index);
		}
	};

#endif // hsp_h
