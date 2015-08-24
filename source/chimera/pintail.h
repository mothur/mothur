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
//This class was created using the algorithms described in the 
// "At Least 1 in 20 16S rRNA Sequence Records Currently Held in the Public Repositories is Estimated To Contain Substantial Anomalies" paper 
//by Kevin E. Ashelford 1, Nadia A. Chuzhanova 3, John C. Fry 1, Antonia J. Jones 2 and Andrew J. Weightman 1.

/***********************************************************/

class Pintail : public Chimera {
	
	public:
		Pintail(string, string, bool, int, string, string, string, int, int, string); //fastafile, templatefile, filter, processors, mask, conservation, quantile, window, increment, outputDir)	
		~Pintail();
		
		int getChimeras(Sequence*);
		Sequence print(ostream&, ostream&);
		
		void setCons(string c)		{ consfile = c;  }
		void setQuantiles(string q) { quanfile = q;  }

	private:
	
		Dist* distcalculator;
		DeCalculator* decalc;
		int iters, window, increment, processors;
		string fastafile, quanfile, consfile;
		
		vector<linePair*> templateLines;
		Sequence* querySeq;
				
		Sequence* bestfit;  //closest match to query in template
		
		vector<float>  obsDistance;  //obsDistance is the vector of observed distances for query 
		vector<float>  expectedDistance;  //expectedDistance is the vector of expected distances for query
		float deviation;  //deviation is the percentage of mismatched pairs over the whole seq between query and its best match.
		vector<int>  windowsForeachQuery;  // windowsForeachQuery is a vector containing the starting spot in query aligned sequence for each window.
										//this is needed so you can move by bases and not just spots in the alignment
										
		int  windowSizes;			//windowSizes = window size of query
		vector<int> windowSizesTemplate;    //windowSizesTemplate[0] = window size of templateSeqs[0]
		
		map<int, int> trimmed;    //trimmed = start and stop of trimmed sequences for query
		map<int, int>::iterator it;
		
		vector<float>  Qav;	//Qav is the vector of average variablility for query
		float  seqCoef;		//seqCoef is the coeff for query
		float DE;			//DE is the deviaation for query
		vector<float> probabilityProfile;
		vector< vector<float> > quantiles;  //quantiles[0] is the vector of deviations with ceiling score of 1, quantiles[1] is the vector of deviations with ceiling score of 2...
		vector< vector<float> > quantilesMembers;  //quantiles[0] is the vector of deviations with ceiling score of 1, quantiles[1] is the vector of deviations with ceiling score of 2...
		set<int>  h;
		string mergedFilterString;
		
		vector< vector<float> > readQuantiles();
		vector<float> readFreq();
		Sequence* findPairs(Sequence*);
			
		void createProcessesQuan();
		int doPrep();
		void printQuanFile(string, string);
		
};

/***********************************************************/

#endif

