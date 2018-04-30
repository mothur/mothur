#ifndef GETLISTCOUNTCOMMAND_H
#define GETLISTCOUNTCOMMAND_H
/*
 *  getlistcountcommand.h
 *  Mothur
 *
 *  Created by westcott on 10/12/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */

#include "command.hpp"
#include "inputdata.h"
#include "listvector.hpp"

/**********************************************************/

class GetListCountCommand : public Command {
	
public:
	GetListCountCommand(string);
	GetListCountCommand();	
	~GetListCountCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "get.otulist";				}
	string getCommandCategory()		{ return "OTU-Based Approaches";	}
	string getCitation() { return "http://www.mothur.org/wiki/Get.otulist"; }
	string getDescription()		{ return "lists each OTU number and the sequence contained in that OTU"; }

	
	string getHelpString();	
    string getOutputPattern(string);	
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
	
private:
	ListVector* list;
	InputData* input;
	
	bool abort, allLines;
	set<string> labels; //holds labels to be used
	string label, listfile, outputDir, sort;
	ofstream out;
	vector<string> outputNames;
		
	void process(ListVector*);
};
/**********************************************************/

#endif
