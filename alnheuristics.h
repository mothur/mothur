#ifndef alnheuristics_h
#define alnheuristics_h

struct AlnParams;

struct AlnHeuristics
	{
	unsigned BandRadius;
	unsigned HSPFinderWordLength;
	float SeedT;

	float XDropG;			//  GappedBlast default
	float XDropU;			//  UngappedBlast default
	float XDropUG;			//  UngappedBlast called by GappedBlast

	unsigned MinGlobalHSPLength;

	AlnHeuristics();
	void InitFromCmdLine(const AlnParams &AP);
	void InitGlobalFull();

	bool IsGlobalFull() const
		{
		return MinGlobalHSPLength == 0 && BandRadius == 0;
		}

	};

#endif // alnheuristics_h
