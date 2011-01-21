#ifndef REMOVEOTUSCOMMAND_H
#define REMOVEOTUSCOMMAND_H

/*
 *  removeotuscommand.h
 *  Mothur
 *
 *  Created by westcott on 11/12/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "command.hpp"
#include "groupmap.h"
#include "listvector.hpp"

class RemoveOtusCommand : public Command {
	
public:
	
	RemoveOtusCommand(string);	
	RemoveOtusCommand();
	~RemoveOtusCommand(){}
	vector<string> getRequiredParameters();
	vector<string> getValidParameters();
	vector<string> getRequiredFiles();
	map<string, vector<string> > getOutputFiles() { return outputTypes; }
	int execute();
	void help();	
	
private:
	string accnosfile, groupfile, listfile, outputDir, groups, label;
	bool abort;
	vector<string> outputNames, Groups;
	map<string, vector<string> > outputTypes;
	GroupMap* groupMap;
	
	void readAccnos();
	int readListGroup();
	int processList(ListVector*&, GroupMap*&, ofstream&, ofstream&, bool&);
	
};

#endif




