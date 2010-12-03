#ifndef INDICATORCOMMAND_H
#define INDICATORCOMMAND_H

/*
 *  indicatorcommand.h
 *  Mothur
 *
 *  Created by westcott on 11/12/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "command.hpp"
#include "readtree.h"
#include "treemap.h"
#include "globaldata.hpp"
#include "sharedrabundvector.h"
#include "sharedrabundfloatvector.h"
#include "inputdata.h"

class IndicatorCommand : public Command {
public:
	IndicatorCommand(string);
	IndicatorCommand();
	~IndicatorCommand();
	vector<string> getRequiredParameters();
	vector<string> getValidParameters();
	vector<string> getRequiredFiles();
	map<string, vector<string> > getOutputFiles() { return outputTypes; }
	int execute();
	void help();
	
private:
	GlobalData* globaldata;
	ReadTree* read;
	TreeMap* treeMap;
	string treefile, sharedfile, relabundfile, groups, label, inputFileName, outputDir;
	bool abort;
	vector<string> outputNames, Groups;
	map<string, vector<string> > outputTypes;
	vector<SharedRAbundVector*> lookup;
	vector<SharedRAbundFloatVector*> lookupFloat;
	
	int getShared();
	int getSharedFloat();
	int GetIndicatorSpecies(Tree*&);
	set<string> getDescendantList(Tree*&, int, map<int, set<string> >, map<int, set<int> >&);
	vector<float> getValues(vector< vector<SharedRAbundVector*> >&);
	vector<float> getValues(vector< vector<SharedRAbundFloatVector*> >&);
	map<int, float> getLengthToLeaf(Tree*&);
	
};


#endif

