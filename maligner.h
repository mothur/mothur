#ifndef MALIGNER_H
#define MALIGNER_H
/*
 *  maligner.h
 *  Mothur
 *
 *  Created by westcott on 9/23/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */
 
#include "decalc.h"

/***********************************************************************/
//This class was modeled after the chimeraMaligner written by the Broad Institute
/***********************************************************************/
struct score_struct {
	int prev;
	int score;
	int row;
	int col;
};
/***********************************************************************/
struct trace_struct {
	int col;
	int oldCol;
	int row;
};
/***********************************************************************/
struct results {
	int regionStart;
	int regionEnd;
	string parent;
	float queryToParent;
	float queryToParentLocal;
	float divR;
};

/**********************************************************************/
class Maligner {

	public:
		
		Maligner(vector<Sequence*>, int, int, int, float, int);
		~Maligner() {};
		
		string getResults(Sequence*);
		float getPercentID() {	return percentIdenticalQueryChimera;	}
		vector<results> getOutput()  {	return outputResults;			}
		
				
	private:
		DeCalculator* decalc;
		Sequence* query;
		vector<Sequence*> refSeqs;
		vector<Sequence*> db;
		int numWanted, matchScore, misMatchPenalty, minCoverage;
		float minDivR, percentIdenticalQueryChimera;
		vector<results> outputResults;
		
		vector<Sequence*> minCoverageFilter(vector<Sequence*>);  //removes top matches that do not have minimum coverage with query.
		int computeChimeraPenalty();
		void verticalFilter(vector<Sequence*>);
		
		vector< vector<score_struct> > buildScoreMatrix(int, int);
		void fillScoreMatrix(vector<vector<score_struct> >&, vector<Sequence*>, int);
		vector<score_struct> extractHighestPath(vector<vector<score_struct> >);
		vector<trace_struct> mapTraceRegionsToAlignment(vector<score_struct>, vector<Sequence*>);
		string constructChimericSeq(vector<trace_struct>, vector<Sequence*>);
		float computePercentID(string, string);
		
};

/***********************************************************************/

#endif

