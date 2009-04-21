#ifndef TREEGROUPCOMMAND_H
#define TREEGROUPCOMMAND_H

/*
 *  treegroupscommand.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 4/8/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */
 
#include "command.hpp"
#include "sharedordervector.h"
#include "sharedlistvector.h"
#include "inputdata.h"
#include "groupmap.h"
#include "readotu.h"
#include "validcalculator.h"
#include "tree.h"
#include "treemap.h"
#include "sharedutilities.h"


/* This command create a tree file for each similarity calculator at distance level, using various calculators to find the similiarity between groups. 
	The user can select the lines or labels they wish to use as well as the groups they would like included.
	They can also use as many or as few calculators as they wish. */
	
class GlobalData;

class TreeGroupCommand : public Command {
	
public:
	TreeGroupCommand();	
	~TreeGroupCommand();
	int execute();	
	
private:
	void createTree();
	void printSims();
	
	GlobalData* globaldata;
	SharedUtil* util;
	ReadOTUFile* read;
	TreeMap* tmap;
	Tree* t;
	vector<Calculator*> treeCalculators;
	vector< vector<float> > simMatrix;
	map<int, int> index;  //maps row in simMatrix to vector index in the tree	
	InputData* input;
	ValidCalculators* validCalculator;
	SharedListVector* SharedList;
	SharedOrderVector* order;
	vector<SharedRAbundVector*> lookup;
	string format, outputFile, groupNames;
	int numGroups;

};
	
	
#endif


