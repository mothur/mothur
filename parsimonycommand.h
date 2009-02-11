#ifndef PARSIMONYCOMMAND_H
#define PARSIMONYCOMMAND_H
/*
 *  parsimonycommand.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 1/26/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <cmath>
#include "command.hpp"
#include "parsimony.h"
#include "treemap.h"

using namespace std;

class GlobalData;

class ParsimonyCommand : public Command {
	
	public:
		ParsimonyCommand();	
		~ParsimonyCommand() { delete pars; }
		int execute();	
	
	private:
		GlobalData* globaldata;
		vector<Tree*> T;	   //user trees
		Tree* randT;  //random tree
		TreeMap* tmap;
		Parsimony* pars;
		string parsFile, sumFile, distFile;
		int iters, randomtree, numGroups;
		vector<int> numEachGroup; //vector containing the number of sequences in each group the users wants for random distrib.
		vector<float> userTreeScores; //scores for users trees
		vector<float> UScoreSig;  //tree score signifigance when compared to random trees - percentage of random trees with that score or lower.
		EstOutput userData;			//pscore info for user tree
		EstOutput randomData;		//pscore info for random trees
		map<int, float> validScores;  //contains scores from both user and random
		map<int, float> rscoreFreq;  //pscore, number of random trees with that score.
		map<int, float> uscoreFreq;  //pscore, number of user trees with that score.
		map<int, float> rCumul;  //pscore, cumulative percentage of number of random trees with that score or lower.
		map<int, float> uCumul;  //pscore, cumulative percentage of number of user trees with that score or lower .
		map<int, float>::iterator it;
		map<int, float>::iterator it2;
		
		ofstream out, outSum, outDist;
		
		void printParsimonyFile();  
		void printUSummaryFile();
		void getUserInput();
		
};


#endif
