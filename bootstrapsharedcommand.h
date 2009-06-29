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
#include "readotu.h"
#include "validcalculator.h"
#include "tree.h"
#include "treemap.h"
#include "sharedutilities.h"
#include "consensuscommand.h"
	
class GlobalData;

class BootSharedCommand : public Command {
	
public:
	BootSharedCommand(string);	
	~BootSharedCommand();
	int execute();	
	void help();
	
private:
	void createTree(ostream*, Tree*);
	void printSims();
	void process(SharedOrderVector*);
	
	
	GlobalData* globaldata;
	SharedUtil* util;
	ReadOTUFile* read;
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
	set<int> lines; //hold lines to be used
	set<string> labels; //holds labels to be used
	string outputFile, calc, groups, line, label;
	int numGroups, iters;
	vector<string>  Estimators, Groups; //holds estimators to be used

};
	
	
#endif


