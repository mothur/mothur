#ifndef READOTUCOMMAND_H
#define READOTUCOMMAND_H
/*
 *  readotu.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 1/20/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "command.hpp"
#include "readotu.h"
#include "inputdata.h"
#include "groupmap.h"
#include "sharedcommand.h"

class ReadOtuCommand : public Command {
public:
	ReadOtuCommand(string);
	ReadOtuCommand();
	~ReadOtuCommand() {}
	
	vector<string> setParameters() {  return outputNames; } //dummy doesn't really do anything
	string getCommandName()			{ return "read.otu";	}
	string getCommandCategory()		{ return "Hidden";	}
	string getHelpString() { return "This command is no longer available. You can provide your files directly to the downstream commands like collect.shared."; }	
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	InputData* input;
	Command* shared;
	GroupMap* groupMap;
	string filename, listfile, orderfile, sharedfile, label, groupfile, sabundfile, rabundfile, format, groups, outputDir, ordergroupfile, relAbundfile;
	vector<string> Groups, outputNames;
	map<string, vector<string> > outputTypes;

	bool abort, allLines;
	set<string> labels; //holds labels to be used

};

#endif
