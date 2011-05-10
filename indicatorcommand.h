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
#include "sharedrabundvector.h"
#include "sharedrabundfloatvector.h"
#include "inputdata.h"

class IndicatorCommand : public Command {
public:
	IndicatorCommand(string);
	IndicatorCommand();
	~IndicatorCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "indicator";				}
	string getCommandCategory()		{ return "Hypothesis Testing";		}
	string getHelpString();	
	string getCitation() { return "http://www.mothur.org/wiki/Indicator"; }
	
	int execute();
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	ReadTree* read;
	TreeMap* treeMap;
	GroupMap* designMap;
	string treefile, sharedfile, relabundfile, groups, label, inputFileName, outputDir, designfile;
	bool abort;
	vector<string> outputNames, Groups;
	vector<SharedRAbundVector*> lookup;
	vector<SharedRAbundFloatVector*> lookupFloat;
	
	int getShared();
	int getSharedFloat();
	int GetIndicatorSpecies(Tree*&);
	set<string> getDescendantList(Tree*&, int, map<int, set<string> >, map<int, set<int> >&);
	vector<float> getValues(vector< vector<SharedRAbundVector*> >&);
	vector<float> getValues(vector< vector<SharedRAbundFloatVector*> >&);
	map<int, float> getDistToRoot(Tree*&);
	
};


#endif

