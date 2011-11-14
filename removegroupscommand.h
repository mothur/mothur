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
	
	vector<string> setParameters();
	string getCommandName()			{ return "remove.groups";			}
	string getCommandCategory()		{ return "OTU-Based Approaches";	}
	string getHelpString();	
	string getCitation() { return "http://www.mothur.org/wiki/Remove.groups"; }
	string getDescription()		{ return "removes sequences from a list, fasta, name, group or taxonomy file from a given group or set of groups"; }
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
	
private:
	set<string> names;
	string accnosfile, fastafile, namefile, groupfile, listfile, taxfile, outputDir, groups, sharedfile;
	bool abort;
	vector<string> outputNames, Groups;
	GroupMap* groupMap;
	map<string, string> uniqueToRedundant; //if a namefile is given and the first column name is not selected
	//then the other files need to change the unique name in their file to match.
	//only add the names that need to be changed to keep the map search quick
	
	
	int readFasta();
	int readShared();
	int readName();
	int readGroup();
	void readAccnos();
	int readList();
	int readTax();
	int fillNames();
	
};

#endif


