#ifndef BOOTSTRAPSHAREDCOMMAND_H
#define BOOTSTRAPSHAREDCOMMAND_H

/*
 *  bootstrapsharedcommand.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 4/16/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */
 
#include "command.hpp"
#include "sharedordervector.h"
#include "inputdata.h"
#include "groupmap.h"
#include "validcalculator.h"
#include "tree.h"
#include "treemap.h"
#include "sharedutilities.h"
#include "consensuscommand.h"

class BootSharedCommand : public Command {
	
public:
	BootSharedCommand(string);	
	BootSharedCommand();
	~BootSharedCommand();
	
	vector<string> setParameters();
	string getCommandName()			{ return "bootstrap.shared";	}
	string getCommandCategory()		{ return "Hidden";				}
	string getHelpString();	
	string getCitation() { return "no citation"; }

	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	int createTree(ostream*, Tree*);
	void printSims();
	int process(SharedOrderVector*);
	
	SharedUtil* util;
	TreeMap* tmap;
	Tree* t;
	Tree* tempTree;
	ConcensusCommand* consensus;
	vector< vector<Tree*> > trees;  //a vector of trees for each calculator chosen
	vector<Calculator*> treeCalculators;
	vector<ofstream*> out;
	vector< vector<float> > simMatrix;
	map<int, int> index;  //maps row in simMatrix to vector index in the tree	
	InputData* input;
	ValidCalculators* validCalculator;
	SharedOrderVector* order;
	vector<SharedRAbundVector*> lookup;

	bool abort, allLines;
	set<string> labels; //holds labels to be used
	string outputFile, calc, groups, label, outputDir, sharedfile;
	int numGroups, iters;
	vector<string>  Estimators, Groups, outputNames; //holds estimators to be used
};
	
	
#endif


