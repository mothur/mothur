#ifndef LIBSHUFFCOMMAND_H
#define LIBSHUFFCOMMAND_H

/*
 *  libshuffcommand.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 3/9/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "command.hpp"
#include "coverage.h"
#include "fullmatrix.h"

using namespace std;

class GlobalData;

class LibShuffCommand : public Command {
	
	public:
		LibShuffCommand();	
		~LibShuffCommand();
		int execute();	
	
	private:
		vector< vector< vector<float> > > cValues; // vector<vector of coverage scores, one for each comparison.> -one for each distance level.
		vector< vector<float> > deltaValues; // vector< vector of delta scores, one for each comparison.> -one at each distance
		vector<float> sumDelta; //sum of delta scores, one for each comparison.
		vector<float> sumDeltaSig; //number of random  matrixes with that delta value or ??
		vector< vector<float> > rsumDelta; //vector< vector<sumdelta scores for a given comparison> >
		vector<string> groupComb;
		vector<float> dist;
		
		
		void setGroups();
		int findIndex(float, int);
		void printCoverageFile();
		void printSummaryFile();

		GlobalData* globaldata;
		Coverage* coverage;
		FullMatrix* matrix;
		float cutOff, step;
		int numGroups, numComp, iters;
		string coverageFile, summaryFile, form;
		ofstream out, outSum;
				
		
};

#endif
