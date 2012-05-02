#ifndef REFCHIMERATEST
#define REFCHIMERATEST

/*
 *  refchimeratest.h
 *  Mothur
 *
 *  Created by Pat Schloss on 1/31/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "sequence.hpp"


class RefChimeraTest {
	
public:
	RefChimeraTest(vector<Sequence>&);
	int printHeader(ofstream&);
	int analyzeQuery(string, string, ofstream&);
	int getClosestRefIndex();
private:
	int getMismatches(string&, vector<vector<int> >&, vector<vector<int> >&, int&);
	int getChimera(vector<vector<int> >&, vector<vector<int> >&, int&, int&, int&, vector<int>&, vector<int>&, vector<int>&, vector<int>&);
	int getTrimera(vector<vector<int> >&, vector<vector<int> >&, int&, int&, int&, int&, int&, vector<int>&, vector<int>&, vector<int>&, vector<int>&);
	string stitchBimera(int, int, int);
	string stitchTrimera(int, int, int, int, int);
	double calcDistToChimera(string&, string&);

	vector<string> referenceSeqs;
	vector<string> referenceNames;
	int numRefSeqs;
	int alignLength;
	int bestMatch;
	//ofstream chimeraReportFile;
	
	MothurOut* m;
};

#endif
