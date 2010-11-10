#ifndef REMOVEGROUPSCOMMAND_H
#define REMOVEGROUPSCOMMAND_H

/*
 *  removegroupscommand.h
 *  Mothur
 *
 *  Created by westcott on 11/10/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */


#include "command.hpp"
#include "groupmap.h"

class RemoveGroupsCommand : public Command {
	
public:
	
	RemoveGroupsCommand(string);	
	RemoveGroupsCommand();
	~RemoveGroupsCommand(){}
	vector<string> getRequiredParameters();
	vector<string> getValidParameters();
	vector<string> getRequiredFiles();
	map<string, vector<string> > getOutputFiles() { return outputTypes; }
	int execute();
	void help();	
	
private:
	set<string> names;
	string accnosfile, fastafile, namefile, groupfile, listfile, taxfile, outputDir, groups;
	bool abort;
	vector<string> outputNames, Groups;
	map<string, vector<string> > outputTypes;
	GroupMap* groupMap;
	
	int readFasta();
	int readName();
	int readGroup();
	void readAccnos();
	int readList();
	int readTax();
	int fillNames();
	
};

#endif


