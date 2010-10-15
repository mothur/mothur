#ifndef PARSELISTCOMMAND_H
#define PARSELISTCOMMAND_H
/*
 *  parselistcommand.h
 *  Mothur
 *
 *  Created by westcott on 2/24/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "command.hpp"
#include "groupmap.h"
#include "inputdata.h"
#include "listvector.hpp"

/***************************************************************************************/

class ParseListCommand : public Command {
	
public:
	ParseListCommand(string);
	ParseListCommand();	
	~ParseListCommand();
	vector<string> getRequiredParameters();
	vector<string> getValidParameters();
	vector<string> getRequiredFiles();
	map<string, vector<string> > getOutputFiles() { return outputTypes; }
	int execute();	
	void help();
	
private:
	int parse(ListVector*);
		
	ListVector* list;
	GroupMap* groupMap;
	InputData* input;
	
	ofstream out;
	string outputDir, listfile, groupfile, label;
	set<string> labels;
	bool abort, allLines;
	vector<string> outputNames;
	map<string, vector<string> > outputTypes;
	
	map<string, ofstream*> filehandles;
	map<string, ofstream*>::iterator it3;

};

/***************************************************************************************/

#endif

