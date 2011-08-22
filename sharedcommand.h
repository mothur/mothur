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
	string getCitation() { return "http://www.mothur.org/wiki/Make.shared"; }
	string getDescription()		{ return "make a shared file from a list and group file"; }

	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	void printSharedData(vector<SharedRAbundVector*>);
	int createMisMatchFile();
	int readOrderFile();
	bool isValidGroup(string, vector<string>);
	int eliminateZeroOTUS(vector<SharedRAbundVector*>&);
	int ListGroupSameSeqs();
	
	SharedListVector* SharedList;
	InputData* input;
	GroupMap* groupMap;
	vector<string> Groups, outputNames, order;
	set<string> labels;
	ofstream out;
	string filename, fileroot, outputDir, listfile, groupfile, ordergroupfile;
	bool firsttime, pickedGroups, abort, allLines;
	map<string, ofstream*> filehandles;
	map<string, ofstream*>::iterator it3;

};

#endif
