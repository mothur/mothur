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
		string weightedFile, sumFile;
		int iters, numGroups, numComp;
		EstOutput userData;			//weighted score info for user tree
		EstOutput randomData;		//weighted score info for random trees
		vector< vector<float> > validScores;  //vector<contains scores from both user and random> each group comb has an entry
		vector< vector<float> > rScores;  //vector<weighted scores for random trees.> each group comb has an entry
		vector< vector<float> > uScores;  //vector<weighted scores for user trees.> each group comb has an entry
								
		ofstream outSum, out;
		
		void printWSummaryFile();
	//	void printWeightedFile();  
		void removeValidScoresDuplicates();
		int findIndex(float);
		void setGroups(); 
};



#endif
