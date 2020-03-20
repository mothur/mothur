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
#include "treereader.h"
#include "parsimony.h"
#include "counttable.h"
#include "fileoutput.h"
#include "readtree.h"


class ParsimonyCommand : public Command {

public:
	ParsimonyCommand(string);	
	~ParsimonyCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "parsimony";				}
	string getCommandCategory()		{ return "Hypothesis Testing";		}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "Slatkin M, Maddison WP (1989). A cladistic measure of gene flow inferred from the phylogenies of alleles. Genetics 123: 603-13. \nSlatkin M, Maddison WP (1990). Detecting isolation by distance using phylogenies of genes. Genetics 126: 249-60. \nMartin AP (2002). Phylogenetic approaches for describing and comparing the diversity of microbial communities. Appl Environ Microbiol 68: 3673-82. \nSchloss PD, Handelsman J (2006). Introducing TreeClimber, a test to compare microbial community structure. Appl Environ Microbiol 72: 2379-84.\nhttp://www.mothur.org/wiki/Parsimony"; }
	string getDescription()		{ return "generic test that describes whether two or more communities have the same structure"; }

	int execute();
	void help() { m->mothurOut(getHelpString()); }
	
private:
	FileOutput* output;
	vector<Tree*> T;	   //user trees
	Tree* randT;  //random tree
	Tree* copyUserTree; 
	CountTable* ct; 
	CountTable* savect;
	vector<string> groupComb; // AB. AC, BC...
	string sumFile, randomtree, allGroups, outputDir, treefile, groupfile, namefile, countfile;
	int iters, numGroups, numComp, counter, processors, numUniquesInName;
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
	
	ofstream outSum;
	bool abort;
	string groups, itersString;
	vector<string> Groups, outputNames; //holds groups to be used
	map<string, string> nameMap;
    vector<string> Treenames; 
	
	void printParsimonyFile();  
	int printUSummaryFile();
	void getUserInput();
	int readNamesFile();
	
};


#endif
