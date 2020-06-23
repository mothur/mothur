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

//**********************************************************************************************************************
struct tidy {
    string otu;
    string group;
    int abund;
    
    tidy() : group(""), otu(""), abund(0) {}
    tidy(string o, string g, int a) : otu(o), group(g), abund(a) {}
};
//**********************************************************************************************************************


/* The shared() command:
	The shared command can only be executed after a successful read.shared command.  
	The shared command parses a .list file and separates it into groups.  
	It outputs a .shared file containing the OTU information for each group.  
	There are no shared command parameters.  The shared command should be in the following format: shared(). */


class SharedCommand : public Command {
	
public:
	SharedCommand(string);	
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
	void printSharedData(SharedRAbundVectors*&, ofstream&, bool&);
	int readOrderFile();
	bool isValidGroup(string, vector<string>);
	int ListGroupSameSeqs(vector<string>&, SharedListVector*);
    int createSharedFromListGroup();
    int createSharedFromBiom();
    int createSharedFromCount();
    void convertSharedFormat();
    string findFormat();
    string getTag(string&);
    vector<string> readRows(string, int&);
    int getDims(string, int&, int&);
    SharedRAbundVectors* readData(string, string, string, vector<string>&, int);
	
	vector<string> Groups, outputNames, order;
	set<string> labels;
	string fileroot,  listfile, groupfile, biomfile, ordergroupfile, countfile, sharedfile;
	bool firsttime, pickedGroups, abort, allLines, keepZeroes;

};

#endif
