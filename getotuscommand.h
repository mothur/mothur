#ifndef GETOTUSCOMMAND_H
#define GETOTUSCOMMAND_H

/*
 *  getotuscommand.h
 *  Mothur
 *
 *  Created by westcott on 11/10/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */



#include "command.hpp"
#include "groupmap.h"
#include "listvector.hpp"

class GetOtusCommand : public Command {
	
public:
	
	GetOtusCommand(string);	
	GetOtusCommand();
	~GetOtusCommand(){}
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



