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
#include "chimera.h"
#include "database.hpp"

/***********************************************************************/
//This class was modeled after the chimeraMaligner written by the Broad Institute
/**********************************************************************/
class Maligner {

	public:
		
		Maligner(vector<Sequence*>, int, int, int, float, int, int, string, Database*, Database*);
		~Maligner() {};
		
		string getResults(Sequence*, DeCalculator*);
		float getPercentID() {	return percentIdenticalQueryChimera;	}
		vector<results> getOutput()  {	return outputResults;			}
		
				
	private:
		Sequence* query;
		vector<Sequence*> refSeqs;
		vector<Sequence*> db;
		int numWanted, matchScore, misMatchPenalty, minCoverage, minSimilarity;
		string searchMethod;
		float minDivR, percentIdenticalQueryChimera;
		vector<results> outputResults;
		vector<int> indexes;  //stores index into template seqs of the refSeqs, so we can return the whole sequence rather than the trimmed and filtered one
		map<int, int> spotMap;
		Database* databaseLeft;
		Database* databaseRight;
		
		vector<Sequence*> minCoverageFilter(vector<Sequence*>);  //removes top matches that do not have minimum coverage with query.
		int computeChimeraPenalty();
		void verticalFilter(vector<Sequence*>);
		
		vector< vector<score_struct> > buildScoreMatrix(int, int);
		void fillScoreMatrix(vector<vector<score_struct> >&, vector<Sequence*>, int);
		vector<score_struct> extractHighestPath(vector<vector<score_struct> >);
		vector<trace_struct> mapTraceRegionsToAlignment(vector<score_struct>, vector<Sequence*>);
		string constructChimericSeq(vector<trace_struct>, vector<Sequence*>);
		float computePercentID(string, string);
		string chimeraMaligner(int, DeCalculator*);
		vector<Sequence*> getBlastSeqs(Sequence*, int);
		vector<Sequence*> getKmerSeqs(Sequence*, int);
		MothurOut* m;
		
};

/***********************************************************************/

#endif

