#include "myutils.h"
#include "chime.h"

void WriteChimeFileHdr(FILE *f)
	{
	if (f == 0)
		return;

	fprintf(f,
		"\tQuery"		// 1
		"\tA"			// 2
		"\tB"			// 3
		"\tIdQM"		// 4
		"\tIdQA"		// 5
		"\tIdQB"		// 6
		"\tIdAB"		// 7
		"\tIdQT"		// 8
		"\tLY"			// 9
		"\tLN"			// 10
		"\tLA"			// 11
		"\tRY"			// 12
		"\tRN"			// 13
		"\tRA"			// 14
		"\tDiv"			// 15
		"\tY"			// 16
		"\n"
		);
	}

void WriteChimeHit(FILE *f, const ChimeHit2 &Hit)
	{
	if (f == 0)
		return;

	if (Hit.Div <= 0.0)
		{
		fprintf(f, "0.0000");		// 0

		fprintf(f,
		  "\t%s", Hit.QLabel.c_str());	// 1

		fprintf(f,
		  "\t*"						// 2
		  "\t*"						// 3
		  "\t*"						// 4
		  "\t*"						// 5
		  "\t*"						// 6
		  "\t*"						// 7
		  "\t*"						// 8
		  "\t*"						// 9
		  "\t*"						// 10
		  "\t*"						// 11
		  "\t*"						// 12
		  "\t*"						// 13
		  "\t*"						// 14
		  "\t*"						// 15
		  "\tN"						// 16
		  "\n"
		  );
		return;
		}

	fprintf(f, "%.4f", Hit.Score);		// 0

	fputc('\t', f);
	fputs(Hit.QLabel.c_str(), f);		// 1

	fputc('\t', f);
	fputs(Hit.ALabel.c_str(), f);		// 2

	fputc('\t', f);
	fputs(Hit.BLabel.c_str(), f);		// 3

	fprintf(f, "\t%.1f", Hit.PctIdQM);	// 4
	fprintf(f, "\t%.1f", Hit.PctIdQA);	// 5
	fprintf(f, "\t%.1f", Hit.PctIdQB);	// 6
	fprintf(f, "\t%.1f", Hit.PctIdAB);	// 7
	fprintf(f, "\t%.1f", Hit.PctIdQT);	// 8

	fprintf(f, "\t%u", Hit.CS_LY);		// 9
	fprintf(f, "\t%u", Hit.CS_LN);		// 10
	fprintf(f, "\t%u", Hit.CS_LA);		// 11

	fprintf(f, "\t%u", Hit.CS_RY);		// 12
	fprintf(f, "\t%u", Hit.CS_RN);		// 13
	fprintf(f, "\t%u", Hit.CS_RA);		// 14

	fprintf(f, "\t%.2f", Hit.Div);		// 15

	fprintf(f, "\t%c", yon(Hit.Accept())); // 16
	fputc('\n', f);
	}

unsigned GetUngappedLength(const byte *Seq, unsigned L)
	{
	unsigned UL = 0;
	for (unsigned i = 0; i < L; ++i)
		if (!isgap(Seq[i]))
			++UL;
	return UL;
	}

void WriteChimeHitX(FILE *f, const ChimeHit2 &Hit)
	{
	if (f == 0)
		return;

	if (Hit.Div <= 0.0)
		return;

	const string &Q3 = Hit.Q3;
	const string &A3 = Hit.A3;
	const string &B3 = Hit.B3;

	const byte *Q3Seq = (const byte *) Q3.c_str();
	const byte *A3Seq = (const byte *) A3.c_str();
	const byte *B3Seq = (const byte *) B3.c_str();

// Aligned
	unsigned ColCount = SIZE(Q3);
	asserta(SIZE(A3) == ColCount && SIZE(B3) == ColCount);

	unsigned LQ = GetUngappedLength(Q3Seq, ColCount);
	unsigned LA = GetUngappedLength(A3Seq, ColCount);
	unsigned LB = GetUngappedLength(B3Seq, ColCount);

	fprintf(f, "\n");
	fprintf(f, "------------------------------------------------------------------------\n");
	fprintf(f, "Query   (%5u nt) %s\n", LQ, Hit.QLabel.c_str());
	fprintf(f, "ParentA (%5u nt) %s\n", LA, Hit.ALabel.c_str());
	fprintf(f, "ParentB (%5u nt) %s\n", LB, Hit.BLabel.c_str());

// Strip terminal gaps in query
	unsigned FromCol = UINT_MAX;
	unsigned ToCol = UINT_MAX;
	for (unsigned Col = 0; Col < ColCount; ++Col)
		{
		if (!isgap(Q3Seq[Col]))
			{
			if (FromCol == UINT_MAX)
				FromCol = Col;
			ToCol = Col;
			}
		}

	unsigned QPos = 0;
	unsigned APos = 0;
	unsigned BPos = 0;
	for (unsigned Col = 0; Col < FromCol; ++Col)
		{
		if (!isgap(A3Seq[Col]))
			++APos;
		if (!isgap(B3Seq[Col]))
			++BPos;
		}

	unsigned Range = ToCol - FromCol + 1;
	unsigned RowCount = (Range + 79)/80;
	unsigned RowFromCol = FromCol;
	for (unsigned RowIndex = 0; RowIndex < RowCount; ++RowIndex)
		{
		fprintf(f, "\n");
		unsigned RowToCol = RowFromCol + 79;
		if (RowToCol > ToCol)
			RowToCol = ToCol;

	// A row
		fprintf(f, "A %5u ", APos + 1);
		for (unsigned Col = RowFromCol; Col <= RowToCol; ++Col)
			{
			char q = Q3Seq[Col];
			char a = A3Seq[Col];
			if (a != q)
				a = tolower(a);
			fprintf(f, "%c", a);
			if (!isgap(a))
				++APos;
			}
		fprintf(f, " %u\n", APos);

	// Q row
		fprintf(f, "Q %5u ", QPos + 1);
		for (unsigned Col = RowFromCol; Col <= RowToCol; ++Col)
			{
			char q = Q3Seq[Col];
			fprintf(f, "%c", q);
			if (!isgap(q))
				++QPos;
			}
		fprintf(f, " %u\n", QPos);

	// B row
		fprintf(f, "B %5u ", BPos + 1);
		for (unsigned Col = RowFromCol; Col <= RowToCol; ++Col)
			{
			char q = Q3Seq[Col];
			char b = B3Seq[Col];
			if (b != q)
				b = tolower(b);
			fprintf(f, "%c", b);
			if (!isgap(b))
				++BPos;
			}
		fprintf(f, " %u\n", BPos);

	// Diffs
		fprintf(f, "Diffs   ");
		for (unsigned Col = RowFromCol; Col <= RowToCol; ++Col)
			{
			char q = Q3Seq[Col];
			char a = A3Seq[Col];
			char b = B3Seq[Col];

			char c = ' ';
			if (isgap(q) || isgap(a) || isgap(b))
				c = ' ';
			else if (Col < Hit.ColXLo)
				{
				if (q == a && q == b)
					c = ' ';
				else if (q == a && q != b)
					c = 'A';
				else if (q == b && q != a)
					c = 'b';
				else if (a == b && q != a)
					c = 'N';
				else
					c = '?';
				}
			else if (Col > Hit.ColXHi)
				{
				if (q == a && q == b)
					c = ' ';
				else if (q == b && q != a)
					c = 'B';
				else if (q == a && q != b)
					c = 'a';
				else if (a == b && q != a)
					c = 'N';
				else
					c = '?';
				}

			fprintf(f, "%c", c);
			}
		fprintf(f, "\n");

	// SNPs
		fprintf(f, "Votes   ");
		for (unsigned Col = RowFromCol; Col <= RowToCol; ++Col)
			{
			char q = Q3Seq[Col];
			char a = A3Seq[Col];
			char b = B3Seq[Col];

			bool PrevGap = Col > 0 && (isgap(Q3Seq[Col-1]) || isgap(A3Seq[Col-1]) || isgap(B3Seq[Col-1]));
			bool NextGap = Col+1 < ColCount && (isgap(Q3Seq[Col+1]) || isgap(A3Seq[Col+1]) || isgap(B3Seq[Col+1]));

			char c = ' ';
			if (isgap(q) || isgap(a) || isgap(b) || PrevGap || NextGap)
				c = ' ';
			else if (Col < Hit.ColXLo)
				{
				if (q == a && q == b)
					c = ' ';
				else if (q == a && q != b)
					c = '+';
				else if (q == b && q != a)
					c = '!';
				else
					c = '0';
				}
			else if (Col > Hit.ColXHi)
				{
				if (q == a && q == b)
					c = ' ';
				else if (q == b && q != a)
					c = '+';
				else if (q == a && q != b)
					c = '!';
				else
					c = '0';
				}

			fprintf(f, "%c", c);
			}
		fprintf(f, "\n");

	// LR row
		fprintf(f, "Model   ");
		for (unsigned Col = RowFromCol; Col <= RowToCol; ++Col)
			{
			if (Col < Hit.ColXLo)
				fprintf(f, "A");
			else if (Col >= Hit.ColXLo && Col <= Hit.ColXHi)
				fprintf(f, "x");
			else
				fprintf(f, "B");
			}

		fprintf(f, "\n");

		RowFromCol += 80;
		}
	fprintf(f, "\n");

	double PctIdBestP = max(Hit.PctIdQA, Hit.PctIdQB);
	double Div = (Hit.PctIdQM - PctIdBestP)*100.0/PctIdBestP;

	unsigned LTot = Hit.CS_LY + Hit.CS_LN + Hit.CS_LA;
	unsigned RTot = Hit.CS_RY + Hit.CS_RN + Hit.CS_RA;

	double PctL = Pct(Hit.CS_LY, LTot);
	double PctR = Pct(Hit.CS_RY, RTot);

	fprintf(f,
	  "Ids.  QA %.1f%%, QB %.1f%%, AB %.1f%%, QModel %.1f%%, Div. %+.1f%%\n",
	  Hit.PctIdQA,
	  Hit.PctIdQB,
	  Hit.PctIdAB,
	  Hit.PctIdQM,
	  Div);

	fprintf(f,
	  "Diffs Left %u: N %u, A %u, Y %u (%.1f%%); Right %u: N %u, A %u, Y %u (%.1f%%), Score %.4f\n",
	  LTot, Hit.CS_LN, Hit.CS_LA, Hit.CS_LY, PctL,
	  RTot, Hit.CS_RN, Hit.CS_RA, Hit.CS_RY, PctR,
	  Hit.Score);
	}
