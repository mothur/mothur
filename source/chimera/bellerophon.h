#ifndef BELLEROPHON_H
#define BELLEROPHON_H

/*
 *  bellerophon.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 7/9/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */


#include "mothurchimera.h"
#include "sparsematrix.hpp"
#include "sequence.hpp"
#include "calculator.h"

typedef list<PCell>::iterator MatData;
typedef map<int, float> SeqMap;  //maps sequence to all distance for that seqeunce

/***********************************************************/

class Bellerophon : public MothurChimera {
	
	public:
		Bellerophon(string, bool, bool, int, int, string);	//fastafile, filter, correction, window, increment, outputDir);
		~Bellerophon() { delete distCalculator; for (int i = 0; i < seqs.size(); i++) { delete seqs[i];  }  seqs.clear(); }
		
		int getChimeras();
		int print(ostream&, ostream&, string);
		
	private:
		DistCalc* distCalculator;
		vector<Sequence*> seqs;
		vector< vector<Preference> > pref; //pref[0] = preference scores for all seqs in window 0.
		string fastafile;
		int iters, count, window, increment, numSeqs; //iters = number of windows
		bool correction;
		
		int generatePreferences(vector<SeqMap>, vector<SeqMap>, int);
		int createSparseMatrix(int, int, SparseMatrix*, vector<Sequence>);
		vector<Preference> getBestPref();
		int driverChimeras(vector<int>);
        int writePrefs(string);
		int readPrefs(string);
		vector<string> getBestWindow();
};

/***********************************************************/

#endif

