#ifndef UNIFRACWEIGHTEDCOMMAND_H
#define UNIFRACWEIGHTEDCOMMAND_H

/*
 *  unifracweightedcommand.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 2/9/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "command.hpp"
#include "weighted.h"
#include "treemap.h"
#include "progress.hpp"
#include "sharedutilities.h"
#include "fileoutput.h"

using namespace std;

class GlobalData;

class UnifracWeightedCommand : public Command {
	
	public:
		UnifracWeightedCommand();	
		~UnifracWeightedCommand() { delete weighted; delete util; }
		int execute();	
	
	private:
		GlobalData* globaldata;
		SharedUtil* util;
		FileOutput* output;
		vector<Tree*> T;	   //user trees
		vector<double> utreeScores;  //user tree unweighted scores
		vector<double> WScoreSig;  //tree weighted score signifigance when compared to random trees - percentage of random trees with that score or lower.
		vector<string> groupComb; // AB. AC, BC...
		Tree* randT;  //random tree
		TreeMap* tmap;
		Weighted* weighted;
		string weightedFile, weightedFileout, sumFile;
		int iters, numGroups, numComp, counter;
		EstOutput userData;			//weighted score info for user tree
		EstOutput randomData;		//weighted score info for random trees
		vector< vector<double> > rScores;  //vector<weighted scores for random trees.> each group comb has an entry
		vector< vector<double> > uScores;  //vector<weighted scores for user trees.> each group comb has an entry
		vector< map<float, float> > rScoreFreq;  //map <weighted score, number of random trees with that score.> -vector entry for each combination.
		vector< map<float, float> > rCumul;  //map <weighted score, cumulative percentage of number of random trees with that score or higher.> -vector entry for each c								
		map<float, float>  validScores;  //map contains scores from random
		map<float, float>::iterator it2;
		map<float, float>::iterator it;
		
		ofstream outSum, out;
		ifstream inFile;
		
		void printWSummaryFile();
		void printWeightedFile();  
		//void removeValidScoresDuplicates();
		int findIndex(float, int);
		void calculateFreqsCumuls();
		//void initFile(string);
		//void output(vector<double>);
		//void resetFile();

};



#endif
