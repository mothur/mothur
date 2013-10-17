#ifndef SHAREDCOMMAND_H
#define SHAREDCOMMAND_H
/*
 *  sharedcommand.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "command.hpp"
#include "sharedlistvector.h"
#include "inputdata.h"

/* The shared() command:
	The shared command can only be executed after a successful read.shared command.  
	The shared command parses a .list file and separates it into groups.  
	It outputs a .shared file containing the OTU information for each group.  
	There are no shared command parameters.  The shared command should be in the following format: shared(). */


class SharedCommand : public Command {
	
public:
	SharedCommand(string);	
	SharedCommand();
	~SharedCommand();
	
	vector<string> setParameters();
	string getCommandName()			{ return "make.shared";				}
	string getCommandCategory()		{ return "OTU-Based Approaches";	}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "http://www.mothur.org/wiki/Make.shared"; }
	string getDescription()		{ return "make a shared file from a list and group file"; }

	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	void printSharedData(vector<SharedRAbundVector*>, ofstream&);
	int readOrderFile();
	bool isValidGroup(string, vector<string>);
	int eliminateZeroOTUS(vector<SharedRAbundVector*>&);
	int ListGroupSameSeqs(vector<string>&, SharedListVector*);
    int createSharedFromListGroup();
    int createSharedFromBiom();
    string getTag(string&);
    vector<string> readRows(string, int&); 
    int getDims(string, int&, int&);
    vector<SharedRAbundVector*> readData(string, string, string, vector<string>&, int);
	
	vector<string> Groups, outputNames, order;
	set<string> labels;
	string fileroot, outputDir, listfile, groupfile, biomfile, ordergroupfile, countfile;
	bool firsttime, pickedGroups, abort, allLines;
	map<string, ofstream*> filehandles;
	map<string, ofstream*>::iterator it3;

};

#endif
