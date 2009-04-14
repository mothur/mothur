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

#include "command.hpp"
#include "parsimony.h"
#include "treemap.h"
#include "progress.hpp"
#include "sharedutilities.h"

using namespace std;

class GlobalData;

class ParsimonyCommand : public Command {
	
	public:
		ParsimonyCommand();	
		~ParsimonyCommand() { delete pars; delete util; }
		int execute();	
	
	private:
		GlobalData* globaldata;
		SharedUtil* util;
		vector<Tree*> T;	   //user trees
		Tree* randT;  //random tree
		Tree* copyUserTree; 
		TreeMap* tmap; 
		TreeMap* savetmap;
		Parsimony* pars;
		vector<string> groupComb; // AB. AC, BC...
		string parsFile, parsFileout, sumFile, randomtree, allGroups;
		int iters, numGroups, numComp, counter;
		vector<int> numEachGroup; //vector containing the number of sequences in each group the users wants for random distrib.
		vector< vector<float> > userTreeScores; //scores for users trees for each comb.
		vector< vector<float> > UScoreSig;  //tree score signifigance when compared to random trees - percentage of random trees with that score or lower.
		EstOutput userData;			//pscore info for user tree
		EstOutput randomData;		//pscore info for random trees
		map<int, double>  validScores;  //map contains scores from both user and random
		vector< map<int, double> > rscoreFreq;  //map <pscore, number of random trees with that score.> -vector entry for each combination.
		vector< map<int, double> > uscoreFreq;  //map <pscore, number of user trees with that score.> -vector entry for each combination.
		vector< map<int, double> > rCumul;  //map <pscore, cumulative percentage of number of random trees with that score or lower.> -vector entry for each combination.
		vector< map<int, double> > uCumul;  //map <pscore, cumulative percentage of number of user trees with that score or lower .> -vector entry for each combination.
		map<int, double>::iterator it;
		map<int, double>::iterator it2;
		
		ofstream out, outSum;
		ifstream inFile;
		
		void printParsimonyFile();  
		void printUSummaryFile();
		void getUserInput();
		void initFile(string);
		void output(vector<double>);
		void resetFile();

		
};


#endif
