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
#include "decalc.h"

/***********************************************************/
//This class was created using the algorythms described in the 
// "At Least 1 in 20 16S rRNA Sequence Records Currently Held in the Public Repositories is Estimated To Contain Substantial Anomalies" paper 
//by Kevin E. Ashelford 1, Nadia A. Chuzhanova 3, John C. Fry 1, Antonia J. Jones 2 and Andrew J. Weightman 1.

/***********************************************************/

class Pintail : public Chimera {
	
	public:
		Pintail(string, string);	
		~Pintail();
		
		int getChimeras();
		void print(ostream&);
		
		void setCons(string c)		{ consfile = c;  }
		void setQuantiles(string q) { quanfile = q;  }
		
		
	private:
	
		Dist* distcalculator;
		DeCalculator* decalc;
		int iters;
		string fastafile, templateFile, consfile;
		
		
		vector<linePair*> lines;
		vector<linePair*> templateLines;
		vector<Sequence*> querySeqs;
		vector<Sequence*> templateSeqs;
		
		vector<Sequence*> bestfit;  //bestfit[0] matches queryseqs[0]...
		
		vector< vector<float> > obsDistance;  //obsDistance[0] is the vector of observed distances for queryseqs[0]... 
		vector< vector<float> > expectedDistance;  //expectedDistance[0] is the vector of expected distances for queryseqs[0]... 
		vector<float> deviation;  //deviation[0] is the percentage of mismatched pairs over the whole seq between querySeqs[0] and its best match.
		vector< vector<int> > windowsForeachQuery;  // windowsForeachQuery[0] is a vector containing the starting spot in queryseqs[0] aligned sequence for each window.
										//this is needed so you can move by bases and not just spots in the alignment
										
		vector<int> windowSizes;			//windowSizes[0] = window size of querySeqs[0]
		vector<int> windowSizesTemplate;    //windowSizesTemplate[0] = window size of templateSeqs[0]
		
		vector< map<int, int> > trimmed;    //trimmed[0] = start and stop of trimmed sequences for querySeqs[0]
		map<int, int>::iterator it;
		
		vector< vector<float> > Qav;	//Qav[0] is the vector of average variablility for queryseqs[0]... 
		vector<float>  seqCoef;				//seqCoef[0] is the coeff for queryseqs[0]...
		vector<float> DE;					//DE[0] is the deviaation for queryseqs[0]...
		vector<float> probabilityProfile;
		vector< vector<float> > quantiles;  //quantiles[0] is the vector of deviations with ceiling score of 1, quantiles[1] is the vector of deviations with ceiling score of 2...
		vector< vector<quanMember> > quantilesMembers;  //quantiles[0] is the vector of deviations with ceiling score of 1, quantiles[1] is the vector of deviations with ceiling score of 2...
		vector< set<int> > h;
		
		
		vector<float> readFreq();
		vector<Sequence*> findPairs(int, int);
			
		void createProcessesSpots();
		void createProcessesPairs();
		void createProcesses();
		void createProcessesQuan();
		
};

/***********************************************************/

#endif

