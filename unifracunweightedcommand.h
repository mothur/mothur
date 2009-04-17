#ifndef UNIFRACUNWEIGHTEDCOMMAND_H
#define UNIFRACUNWEIGHTEDCOMMAND_H

/*
 *  unifracunweightedcommand.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 2/9/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "command.hpp"
#include "unweighted.h"
#include "treemap.h"
#include "sharedutilities.h"
#include "fileoutput.h"


using namespace std;

class GlobalData;

class UnifracUnweightedCommand : public Command {
	
	public:
		UnifracUnweightedCommand();	
		~UnifracUnweightedCommand() { delete unweighted; delete util; }
		int execute();	
	
	private:
		GlobalData* globaldata;
		SharedUtil* util;
		FileOutput* output;
		vector<Tree*> T;	   //user trees
		TreeMap* tmap;
		Unweighted* unweighted;
		string sumFile, allGroups;
		vector<string> groupComb; // AB. AC, BC...
		int iters, numGroups, numComp, counter;
		EstOutput userData;			//unweighted score info for user tree
		EstOutput randomData;		//unweighted score info for random trees
		vector< vector<float> > utreeScores; //scores for users trees for each comb.
		vector< vector<float> > UWScoreSig;  //tree score signifigance when compared to random trees - percentage of random trees with that score or higher.
		map<float, float>  validScores;  //map contains scores from random
		vector< map<float, float> > rscoreFreq;  //map <unweighted score, number of random trees with that score.> -vector entry for each combination.
		vector< map<float, float> > rCumul;  //map <unweighted score, cumulative percentage of number of random trees with that score or higher.> -vector entry for each combination.
		map<float, float>::iterator it2;
		map<float, float>::iterator it;
		
		ofstream outSum, out;
		ifstream inFile;
		
		void printUWSummaryFile();
		void printUnweightedFile();
		 
		
};



#endif