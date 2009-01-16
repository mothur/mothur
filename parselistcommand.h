#ifndef PARSELISTCOMMAND_H
#define PARSELISTCOMMAND_H
/*
 *  parselistcommand.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include <Carbon/Carbon.h>
#include <iostream>
#include <fstream>
#include <map>
#include "command.hpp"
#include "rabundvector.hpp"
#include "listvector.hpp"
#include "inputdata.h"
#include "groupmap.h"
#include "readmatrix.hpp"


/* The parselist() command:
	The parselist command is similar to the shared command.  
	It parses a list file and separates it into groups.   It outputs a .list file for each group.  
	The parselist command parameter options are listfile and groupfile.  
	The parselist command should be in the following format: parselist(listfile=yourListFile, groupfile=yourGroupFile). 
	The listfile parameter and groupfile paramater are required.  */

class GlobalData;

class ParseListCommand : public Command {
	
public:
	ParseListCommand();	
	~ParseListCommand();
	int execute();	
	
private:
	GlobalData* globaldata;
	GroupMap* groupMap;
	InputData* input;
	ReadMatrix* read;
	map<string, ofstream*> filehandles;
	map<string, ListVector*> groupOfLists;
	ListVector* list;
	map<string, string> listGroups; //maps group name to sequences from that group in a specific OTU
	map<string, string>::iterator it;
	map<string, ListVector*>::iterator it2;
	map<string, ofstream*>::iterator it3;
	void parse(int);
	string fileroot;
};

#endif