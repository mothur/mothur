#ifndef PINTAIL_H
#define PINTAIL_H

/*
 *  pintail.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 7/9/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "chimera.h"
#include "dist.h"

//This class was created using the algorythms described in the 
// "At Least 1 in 20 16S rRNA Sequence Records Currently Held in the Public Repositories is Estimated To Contain Substantial Anomalies" paper 
//by Kevin E. Ashelford 1, Nadia A. Chuzhanova 3, John C. Fry 1, Antonia J. Jones 2 and Andrew J. Weightman 1.

/***********************************************************/

class Pintail : public Chimera {
	
	public:
		Pintail(string, string);	
		~Pintail();
		
		void getChimeras();
		void print(ostream&);
		
		void setCons(string c) { consfile = c;  }
		
		
	private:
	
		struct linePair {
			int start;
			int end;
			linePair(int i, int j) : start(i), end(j) {}
		};

		Dist* distcalculator;
		int iters;
		string fastafile, templateFile, consfile;
		vector<linePair*> lines;
		vector<Sequence*> querySeqs;
		vector<Sequence*> templateSeqs;
		
		vector<Sequence> bestfit;  //bestfit[0] matches queryseqs[0]...
		
		vector< vector<float> > obsDistance;  //obsDistance[0] is the vector of observed distances for queryseqs[0]... 
		vector< vector<float> > expectedDistance;  //expectedDistance[0] is the vector of expected distances for queryseqs[0]... 
		vector<float> deviation;  //deviation[0] is the percentage of mismatched pairs over the whole seq between querySeqs[0] and its best match.
		vector< vector<int> > windows;  // windows[0] is a vector containing the starting spot in queryseqs[0] aligned sequence for each window.
										//this is needed so you can move by bases and not just spots in the alignment
		vector< map<int, int> >  trim;  //trim[0] is the start and end position of trimmed querySeqs[0].  Used to find the variability over each sequence window.
										
		vector<int> windowSizes;    //windowSizes[0] = window size of querySeqs[0]
		
		vector< vector<float> > Qav;	//Qav[0] is the vector of average variablility for queryseqs[0]... 
		vector<float>  seqCoef;				//seqCoef[0] is the coeff for queryseqs[0]...
		vector<float> DE;					//DE[0] is the deviaation for queryseqs[0]...
		vector<float> probabilityProfile;
		
		vector<Sequence*> readSeqs(string);
		void trimSeqs(Sequence*, Sequence&, int);
		vector<float> readFreq();
		vector< vector<float> > findQav(int, int);  
		vector<float> calcFreq(vector<Sequence*>);
		vector<float> getCoef(int, int);
		
		vector<Sequence> findPairs(int, int);
		vector< vector<int> > findWindows(int, int);
		vector< vector<float> > calcObserved(int, int);
		vector< vector<float> > calcExpected(int, int);
		vector<float> calcDE(int, int);
		vector<float> calcDist(int, int);
	
		void createProcessesSpots();
		void createProcesses();
		
		
		
};

/***********************************************************/

#endif

