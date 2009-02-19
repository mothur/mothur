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

using namespace std;

class GlobalData;

class UnifracUnweightedCommand : public Command {
	
	public:
		UnifracUnweightedCommand();	
		~UnifracUnweightedCommand() { delete unweighted; }
		int execute();	
	
	private:
		GlobalData* globaldata;
		vector<Tree*> T;	   //user trees
		vector<float> utreeScores;  //user tree unweighted scores
		vector<float> UWScoreSig;  //tree unweighted score signifigance when compared to random trees - percentage of random trees with that score or lower.
		Tree* randT;  //random tree
		TreeMap* tmap;
		Unweighted* unweighted;
		string sumFile, distFile, unweightedFile;
		int iters;
		EstOutput userData;			//unweighted score info for user tree
		EstOutput randomData;		//unweighted score info for random trees
		map<float, float> validScores;  //contains scores from both user and random
		map<float, float> rscoreFreq;  //unweighted score, number of random trees with that score.
		map<float, float> uscoreFreq;  //unweighted, number of user trees with that score.
		map<float, float> totalrscoreFreq;  //unweighted score, number of random trees with that score.
		map<float, float> rCumul;		//unweighted score, cumulative percentage of number of random trees with that score or higher.
		map<float, float> uCumul;  //unweighted, cumulative percentage of number of user trees with that score or higher .
		map<float, float>::iterator it;
		map<float, float>::iterator it2;
		
		ofstream outSum, outDist, out;
		
		void printUWSummaryFile();
		void printUnweightedFile();
		void saveRandomScores();    
		
};



#endif