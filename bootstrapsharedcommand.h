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
#include "sharedlistvector.h"
#include "inputdata.h"
#include "groupmap.h"
#include "readmatrix.hpp"
#include "validcalculator.h"
#include "tree.h"
#include "treemap.h"
#include "sharedutilities.h"
	
class GlobalData;

class BootSharedCommand : public Command {
	
public:
	BootSharedCommand();	
	~BootSharedCommand();
	int execute();	
	
private:
	void createTree(ostream*);
	void printSims();
	
	GlobalData* globaldata;
	SharedUtil* util;
	ReadMatrix* read;
	TreeMap* tmap;
	Tree* t;
	vector<Calculator*> treeCalculators;
	vector<ofstream*> out;
	vector< vector<float> > simMatrix;
	map<int, int> index;  //maps row in simMatrix to vector index in the tree	
	InputData* input;
	ValidCalculators* validCalculator;
	SharedListVector* SharedList;
	SharedOrderVector* order;
	vector<SharedRAbundVector*> lookup;
	string format, outputFile;
	int numGroups, iters;

};
	
	
#endif


