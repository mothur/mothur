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
#include "readtree.h"

class UnifracWeightedCommand : public Command {
	
	public:
		UnifracWeightedCommand(string);
		UnifracWeightedCommand();
		~UnifracWeightedCommand() {}
	
		vector<string> setParameters();
		string getCommandName()			{ return "unifrac.weighted";		}
		string getCommandCategory()		{ return "Hypothesis Testing";		}
		string getHelpString();	
		string getCitation() { return "Lozupone CA, Hamady M, Kelley ST, Knight R (2007). Quantitative and qualitative beta diversity measures lead to different insights into factors that structure microbial communities. Appl Environ Microbiol 73: 1576-85. \nhttp://www.mothur.org/wiki/Unifrac.weighted"; }
		string getDescription()		{ return "generic tests that describes whether two or more communities have the same structure"; }

		int execute();
		void help() { m->mothurOut(getHelpString()); }
	
	private:
		struct linePair {
			int start;
			int num;
			linePair(int i, int j) : start(i), num(j) {}
		};
		vector<linePair> lines;
        TreeMap* tmap;
		FileOutput* output;
		vector<Tree*> T;	   //user trees
		vector<double> utreeScores;  //user tree unweighted scores
		vector<double> WScoreSig;  //tree weighted score signifigance when compared to random trees - percentage of random trees with that score or lower.
		vector<string> groupComb; // AB. AC, BC...
		string sumFile, outputDir;
		int iters, numGroups, numComp, counter;
		vector< vector<double> > rScores;  //vector<weighted scores for random trees.> each group comb has an entry
		vector< vector<double> > uScores;  //vector<weighted scores for user trees.> each group comb has an entry
		vector< map<float, float> > rScoreFreq;  //map <weighted score, number of random trees with that score.> -vector entry for each combination.
		vector< map<float, float> > rCumul;  //map <weighted score, cumulative percentage of number of random trees with that score or higher.> -vector entry for each c								
		map<float, float>  validScores;  //map contains scores from random
		
		bool abort, phylip, random, includeRoot, subsample, consensus;
		string groups, itersString, outputForm, treefile, groupfile, namefile;
		vector<string> Groups, outputNames; //holds groups to be used
		int processors, numUniquesInName, subsampleSize, subsampleIters;
		ofstream outSum;
		map<string, string> nameMap;
		
		void printWSummaryFile();
		void printWeightedFile();  
		void createPhylipFile();
		//void removeValidScoresDuplicates();
		int findIndex(float, int);
		void calculateFreqsCumuls();
		int createProcesses(Tree*,  vector< vector<string> >,  vector< vector<double> >&);
		int driver(Tree*, vector< vector<string> >, int, int,  vector< vector<double> >&);
        int runRandomCalcs(Tree*, vector<double>);
        vector<Tree*> buildTrees(vector< vector<double> >&, int, TreeMap&);
        int getConsensusTrees(vector< vector<double> >&, int);
        int getAverageSTDMatrices(vector< vector<double> >&, int);
		
};



#endif
