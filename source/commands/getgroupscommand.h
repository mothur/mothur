#ifndef GETGROUPSCOMMAND_H
#define GETGROUPSCOMMAND_H

/*
 *  getgroupscommand.h
 *  Mothur
 *
 *  Created by westcott on 11/10/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */


#include "command.hpp"
#include "groupmap.h"

class GetGroupsCommand : public Command {
    
#ifdef UNIT_TEST
    friend class TestGetGroupsCommand;
#endif
	
public:
	
	GetGroupsCommand(string);	
	GetGroupsCommand();
	~GetGroupsCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "get.groups";				}
	string getCommandCategory()		{ return "OTU-Based Approaches";	}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "http://www.mothur.org/wiki/Get.groups"; }
	string getDescription()		{ return "gets sequences from a list, fasta, name, group, shared, design or taxonomy file from a given group or set of groups"; }

	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
	
private:
	set<string> names;
	map<string, string> uniqueToRedundant; //if a namefile is given and the first column name is not selected
										   //then the other files need to change the unique name in their file to match.
										   //only add the names that need to be changed to keep the map search quick
	string accnosfile, countfile, fastafile, namefile, groupfile, listfile, designfile, taxfile, outputDir, groups, sharedfile, phylipfile, columnfile;
	bool abort;
	vector<string> outputNames, Groups;
	GroupMap* groupMap;
	
	int readFasta();
	int readName();
	int readGroup();
    int readCount();
	int readList();
	int readTax();
	int fillNames();
	int readShared();
    int readDesign();
    int readPhylip();
    int readColumn();
	
};

#endif



