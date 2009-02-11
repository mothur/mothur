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

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <cmath>
#include "command.hpp"
#include "weighted.h"
#include "treemap.h"

using namespace std;

class GlobalData;

class UnifracWeightedCommand : public Command {
	
	public:
		UnifracWeightedCommand();	
		~UnifracWeightedCommand() { delete weighted; }
		int execute();	
	
	private:
		GlobalData* globaldata;
		vector<Tree*> T;	   //user trees
		vector<float> utreeScores;  //user tree unweighted scores
		vector<float> WScoreSig;  //tree weighted score signifigance when compared to random trees - percentage of random trees with that score or lower.
		vector<string> groupComb; // AB. AC, BC...
		Tree* randT;  //random tree
		TreeMap* tmap;
		Weighted* weighted;
		string weightedFile, sumFile, distFile;
		int iters, numGroups, numComp;
		EstOutput userData;			//weighted score info for user tree
		EstOutput randomData;		//weighted score info for random trees
		vector< map<float, float> > validScores;  //vector<contains scores from both user and random> each group comb has an entry
		vector< map<float, float> > rscoreFreq;  //vector<weighted score, number of random trees with that score.> each group comb has an entry
		vector< map<float, float> > uscoreFreq;  //vector<weighted, number of user trees with that score.> each group comb has an entry
		vector< map<float, float> > totalrscoreFreq;  //vector<weighted score, number of random trees with that score.> each group comb has an entry
		vector< map<float, float> > rCumul;  //vector<weighted score, number of random trees with that score.> each group comb has an entry
		vector< map<float, float> > uCumul;  //vector<weighted, cumulative percentage of number of user trees with that score or lower.> each group comb has an entry
		map<float, float>::iterator it;
		map<float, float>::iterator it2;
		
		ofstream outSum, outDist, out;
		
		void printWSummaryFile();
		void printWeightedFile();  
		void saveRandomScores(); 
};



#endif
