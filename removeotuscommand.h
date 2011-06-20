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
	
	vector<string> setParameters();
	string getCommandName()			{ return "remove.otus";				}
	string getCommandCategory()		{ return "OTU-Based Approaches";	}
	string getHelpString();	
	string getCitation() { return "http://www.mothur.org/wiki/Remove.otus"; }
	string getDescription()		{ return "outputs a new list file containing the otus NOT containing sequences from the groups specified"; }

	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	string accnosfile, groupfile, listfile, outputDir, groups, label;
	bool abort;
	vector<string> outputNames, Groups;
	GroupMap* groupMap;
	
	void readAccnos();
	int readListGroup();
	int processList(ListVector*&, GroupMap*&, ofstream&, ofstream&, bool&);
	
};

#endif




