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
    
#ifdef UNIT_TEST
    friend class TestRemoveGroupsCommand;
#endif
	
public:
	
	RemoveGroupsCommand(string);	
	~RemoveGroupsCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "remove.groups";			}
	string getCommandCategory()		{ return "OTU-Based Approaches";	}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "http://www.mothur.org/wiki/Remove.groups"; }
	string getDescription()		{ return "removes sequences from a list, fasta, name, group, shared, design or taxonomy file from a given group or set of groups"; }
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
	
private:
	set<string> names;
	string accnosfile, fastafile, namefile, groupfile, countfile, designfile, listfile, taxfile, outputDir, groups, sharedfile, phylipfile, columnfile, sets;
	bool abort;
	vector<string> outputNames, Groups, Sets;
	GroupMap* groupMap;
	map<string, string> uniqueToRedundant; //if a namefile is given and the first column name is not selected
	//then the other files need to change the unique name in their file to match.
	//only add the names that need to be changed to keep the map search quick

    
    void readFasta();
    void readName();
    void readGroup();
    void readList();
    void readTax();
    void fillNames();
    void readShared();
    void readDesign();
    void readPhylip();
    void readColumn();
    void fillGroupsFromDesign();
	
};

#endif


