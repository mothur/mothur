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


#include "chimera.h"
#include "sparsematrix.hpp"
#include "sequence.hpp"
#include "dist.h"

typedef list<PCell>::iterator MatData;
typedef map<int, float> SeqMap;  //maps sequence to all distance for that seqeunce

/***********************************************************/

class Bellerophon : public Chimera {
	
	public:
		Bellerophon(string, bool, bool, int, int, int, string);	//fastafile, filter, correction, window, increment, processors, outputDir);	
		~Bellerophon() { delete distCalculator; for (int i = 0; i < seqs.size(); i++) { delete seqs[i];  }  seqs.clear(); }
		
		int getChimeras();
		int print(ostream&, ostream&, string);
		
	private:
		struct linePair {
			unsigned long long start;
			int num;
			linePair(unsigned long long i, int j) : start(i), num(j) {}
		};
		
		vector<linePair> lines;
	
		Dist* distCalculator;
		vector<Sequence*> seqs;
		vector< vector<Preference> > pref; //pref[0] = preference scores for all seqs in window 0.
		string fastafile;
		int iters, count, window, increment, numSeqs, processors; //iters = number of windows
		bool correction;
		
		int generatePreferences(vector<SeqMap>, vector<SeqMap>, int);
		int createSparseMatrix(int, int, SparseMatrix*, vector<Sequence>);
		vector<Preference> getBestPref();
		int driverChimeras(vector<int>, linePair);
		int createProcesses(vector<int>);
		int writePrefs(string, linePair);
		int readPrefs(string);
		vector<string> getBestWindow(linePair line);
		int fillPref(int, vector<string>&);
};

/***********************************************************/

#endif

