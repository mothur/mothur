#ifndef GETSABUNDCOMMAND_H
#define GETSABUNDCOMMAND_H

/*
 *  getsabundcommand.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 6/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */


#include "command.hpp"
#include "inputdata.h"
#include "sabundvector.hpp"

class GetSAbundCommand : public Command {
public:
	GetSAbundCommand(string);
	GetSAbundCommand();
	~GetSAbundCommand() {}
	
	vector<string> setParameters();
	string getCommandName()			{ return "get.sabund";				}
	string getCommandCategory()		{ return "OTU-Based Approaches";	}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "http://www.mothur.org/wiki/Get.sabund"; }
	string getDescription()		{ return "creates a sabund file"; }

	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	string filename, format, inputfile, listfile, rabundfile, outputDir;
	ofstream out;
	InputData* input;
	SAbundVector* sabund;
	vector<string> outputNames;

	bool abort, allLines;
	set<string> labels; //holds labels to be used
	string label;

};

#endif
