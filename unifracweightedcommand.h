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
#include "counttable.h"
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
    string getOutputPattern(string);	
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
        CountTable* ct;
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
		string groups, itersString, outputForm, treefile, groupfile, namefile, countfile;
		vector<string> Groups, outputNames; //holds groups to be used
		int processors, subsampleSize, subsampleIters;
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
        vector<Tree*> buildTrees(vector< vector<double> >&, int, CountTable&);
        int getConsensusTrees(vector< vector<double> >&, int);
        int getAverageSTDMatrices(vector< vector<double> >&, int);
		
};

/***********************************************************************/
struct weightedRandomData {
    int start;
	int num;
	MothurOut* m;
    vector< vector<double> > scores;
    vector< vector<string> > namesOfGroupCombos;
    Tree* t;
    CountTable* ct;
    bool includeRoot;
	
	weightedRandomData(){}
	weightedRandomData(MothurOut* mout, int st, int en, vector< vector<string> > ngc, Tree* tree, CountTable* count, bool ir, vector< vector<double> > sc) {
        m = mout;
		start = st;
		num = en;
        namesOfGroupCombos = ngc;
        t = tree;
        ct = count;
        includeRoot = ir;
        scores = sc;
	}
};

/**************************************************************************************************/
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
#else
static DWORD WINAPI MyWeightedRandomThreadFunction(LPVOID lpParam){
	weightedRandomData* pDataArray;
	pDataArray = (weightedRandomData*)lpParam;
	try {
        
        Tree* randT = new Tree(pDataArray->ct);
        
        Weighted weighted(pDataArray->includeRoot);
        
		for (int h = pDataArray->start; h < (pDataArray->start+pDataArray->num); h++) {
            
			if (pDataArray->m->control_pressed) { return 0; }
            
			//initialize weighted score
			string groupA = pDataArray->namesOfGroupCombos[h][0];
			string groupB = pDataArray->namesOfGroupCombos[h][1];
			
			//copy T[i]'s info.
			randT->getCopy(pDataArray->t);
            
			//create a random tree with same topology as T[i], but different labels
			randT->assembleRandomUnifracTree(groupA, groupB);
			
			if (pDataArray->m->control_pressed) { delete randT;  return 0;  }
            
			//get wscore of random tree
			EstOutput randomData = weighted.getValues(randT, groupA, groupB);
            
			if (pDataArray->m->control_pressed) { delete randT;  return 0;  }
            
			//save scores
			pDataArray->scores[h].push_back(randomData[0]);
		}
        
		delete randT;
        
        return 0;
    }
	catch(exception& e) {
		pDataArray->m->errorOut(e, "Weighted", "MyWeightedRandomThreadFunction");
		exit(1);
	}
}
#endif


#endif
