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
#include "parselistcommand.h"

/* The read.otu must be run before you execute a collect.single, rarefaction.single, summary.single, 
collect.shared, rarefaction.shared or summary.shared command. Mothur will generate a .list, .rabund and .sabund 
upon completion of the cluster command or you may use your own. The read.otu command parameter options are 
listfile, rabundfile, sabundfile, groupfile and orderfile. The reaad.otu command can be used in two ways. 
The first is to read a listfile, rabundfile or sabundfile and run the collect.single, rarefaction.single or summary.single. 
For this use the read.otu command should be in the following format: read.otu(listfile=yourListFile, orderfile=yourOrderFile). 
The listfile, rabundfile or sabundfile parameter is required, but you may only use one of them. 
The second way to use the read.otu command is to read a listfile and a groupfile so you can use the collect.shared, 
rarefaction.shared or summary.shared commands. In this case the read.otu command should be in the following format: 
read.otu(listfile=yourListFile, groupfile=yourGroupFile). The listfile parameter and groupfile paramaters are required. 
When using the command the second way read.otu command parses the .list file and separates it into groups. 
It outputs a .shared file containing the OTU information for each group. The read.otu command also outputs a .list file for each group. */

class GlobalData;

class ReadOtuCommand : public Command {
public:
	ReadOtuCommand(string);
	~ReadOtuCommand();
	int execute();
	void help();
	
private:
	GlobalData* globaldata;
	ReadOTUFile* read;
	InputData* input;
	Command* shared;
	Command* parselist;
	GroupMap* groupMap;
	string filename, listfile, orderfile, sharedfile, line, label, groupfile, sabundfile, rabundfile, format;
	OptionParser* parser;
	map<string, string> parameters;
	map<string, string>::iterator it;
	bool abort, allLines;
	set<int> lines; //hold lines to be used
	set<string> labels; //holds labels to be used

};

#endif
