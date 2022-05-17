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
#include "counttable.h"

#include "fileoutput.h"
#include "readtree.h"

class UnifracUnweightedCommand : public Command {
	
	public:
		UnifracUnweightedCommand(string);	
		~UnifracUnweightedCommand() = default;
	
		vector<string> setParameters();
		string getCommandName()			{ return "unifrac.unweighted";		}
		string getCommandCategory()		{ return "Hypothesis Testing";		}
		
	string getHelpString();	
    string getOutputPattern(string);	
		string getCitation() { return "Lozupone C, Knight R (2005). UniFrac: a new phylogenetic method for comparing microbial communities. Appl Environ Microbiol 71: 8228-35. \nhttp://www.mothur.org/wiki/Unifrac.unweighted"; }
		string getDescription()		{ return "generic tests that describes whether two or more communities have the same structure"; }

		int execute();
		void help() { m->mothurOut(getHelpString()); }
	
	
	private:

		string sumFile, allGroups;
		vector<string> groupComb; // AB. AC, BC...
		int iters, numGroups, numComp, counter, processors, subsampleSize, subsampleIters, withReplacement;
		vector< vector<float> > utreeScores; //scores for users trees for each comb.
		vector< vector<float> > UWScoreSig;  //tree score signifigance when compared to random trees - percentage of random trees with that score or higher.
		map<float, float>  validScores;  //map contains scores from random
		vector< map<float, float> > rscoreFreq;  //map <unweighted score, number of random trees with that score.> -vector entry for each combination.
		vector< map<float, float> > rCumul;  //map <unweighted score, cumulative percentage of number of random trees with that score or higher.> -vector entry for each combination.
		
		bool abort, phylip, random, includeRoot, consensus, subsample;
		string groups, itersString,  outputForm, treefile, groupfile, namefile, countfile;
		vector<string> Groups, outputNames; //holds groups to be used

		ofstream outSum, out;
		ifstream inFile;
		
        int runRandomCalcs(Tree*, vector<double>);
		void printUWSummaryFile(int);
		void printUnweightedFile(int);
		void createPhylipFile(int);
        vector<Tree*> buildTrees(vector< vector<double> >&, int, CountTable&);
        int getConsensusTrees(vector< vector<double> >&, int);
        int getAverageSTDMatrices(vector< vector<double> >&, int);
		
};

#endif
