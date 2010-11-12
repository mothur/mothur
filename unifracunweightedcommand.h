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


class GlobalData;

class UnifracUnweightedCommand : public Command {
	
	public:
		UnifracUnweightedCommand(string);	
		UnifracUnweightedCommand();
		~UnifracUnweightedCommand() { globaldata->Groups.clear();  if (abort == false) { delete unweighted; delete util; } }
		vector<string> getRequiredParameters();
		vector<string> getValidParameters();
		vector<string> getRequiredFiles();
		map<string, vector<string> > getOutputFiles() { return outputTypes; }
		int execute();
		void help();	
	
	private:
		GlobalData* globaldata;
		SharedUtil* util;
		FileOutput* output;
		vector<Tree*> T;	   //user trees
		TreeMap* tmap;
		Unweighted* unweighted;
		string sumFile, allGroups;
		vector<string> groupComb; // AB. AC, BC...
		int iters, numGroups, numComp, counter, processors;
		EstOutput userData;			//unweighted score info for user tree
		EstOutput randomData;		//unweighted score info for random trees
		vector< vector<float> > utreeScores; //scores for users trees for each comb.
		vector< vector<float> > UWScoreSig;  //tree score signifigance when compared to random trees - percentage of random trees with that score or higher.
		map<float, float>  validScores;  //map contains scores from random
		vector< map<float, float> > rscoreFreq;  //map <unweighted score, number of random trees with that score.> -vector entry for each combination.
		vector< map<float, float> > rCumul;  //map <unweighted score, cumulative percentage of number of random trees with that score or higher.> -vector entry for each combination.
		
		bool abort, phylip, random;
		string groups, itersString, outputDir, outputForm;
		vector<string> Groups, outputNames; //holds groups to be used
		map<string, vector<string> > outputTypes;

		ofstream outSum, out;
		ifstream inFile;
		
		void printUWSummaryFile(int);
		void printUnweightedFile();
		void createPhylipFile(int);
		 
		
};

#endif
