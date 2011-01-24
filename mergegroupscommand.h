#ifndef MERGEGROUPSCOMMAND_H
#define MERGEGROUPSCOMMAND_H

/*
 *  mergegroupscommand.h
 *  mothur
 *
 *  Created by westcott on 1/24/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "command.hpp"
#include "inputdata.h"
#include "sharedrabundvector.h"

class GlobalData;

class MergeGroupsCommand : public Command {
	
public:
	MergeGroupsCommand(string);
	MergeGroupsCommand();
	~MergeGroupsCommand();
	vector<string> getRequiredParameters();
	vector<string> getValidParameters();
	vector<string> getRequiredFiles();
	map<string, vector<string> > getOutputFiles() { return outputTypes; }
	int execute();
	void help();
	
private:
	GlobalData* globaldata;
	GroupMap* designMap;
	vector<SharedRAbundVector*> lookup;
	map<string, vector<string> > outputTypes;
	
	bool abort, allLines, pickedGroups;
	set<string> labels; //holds labels to be used
	string groups, label, outputDir, inputDir, designfile, sharedfile;
	vector<string> Groups, outputNames;
		
	int process(vector<SharedRAbundVector*>&, ofstream&);
};

#endif

