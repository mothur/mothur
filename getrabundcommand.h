#ifndef GETRABUNDCOMMAND_H
#define GETRABUNDCOMMAND_H

/*
 *  getrabundcommand.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 6/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */


#include "command.hpp"
#include "inputdata.h"
#include "listvector.hpp"


class GetRAbundCommand : public Command {
public:
	GetRAbundCommand(string);
	GetRAbundCommand();
	~GetRAbundCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "get.rabund";				}
	string getCommandCategory()		{ return "OTU-Based Approaches";	}
	string getHelpString();	
	string getCitation() { return "http://www.mothur.org/wiki/Get.rabund"; }
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
	
private:
	
	string filename, listfile, sabundfile, inputfile, format, outputDir;
	ofstream out;
	InputData* input;
	RAbundVector* rabund;
	vector<string> outputNames;

	bool abort, allLines, sorted;
	set<string> labels; //holds labels to be used
	string label;

	
};

#endif

