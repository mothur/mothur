#ifndef CATCHALLCOMMAND_H
#define CATCHALLCOMMAND_H

/*
 *  catchallcommand.h
 *  Mothur
 *
 *  Created by westcott on 5/11/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "command.hpp"
#include "inputdata.h"
#include "readotu.h"
#include "sabundvector.hpp"

/* 
 citation goes here
 */ 

/****************************************************************************/

class CatchAllCommand : public Command {

public:

	CatchAllCommand(string);
	CatchAllCommand();
	~CatchAllCommand() {}
	vector<string> getRequiredParameters();
	vector<string> getValidParameters();
	vector<string> getRequiredFiles();
	map< string, vector<string> > getOutputFiles() { return outputTypes; }
	int execute();
	void help();	
	
private:

	GlobalData* globaldata;
	ReadOTUFile* read;
	InputData* input;
	
	string outputDir, sabundfile, rabundfile, listfile, format, path;
	bool abort, allLines;
	set<string> labels;
	vector<string> outputNames;
	map< string, vector<string> > outputTypes;
	
	string process(SAbundVector*);
	int createSummaryFile(string, string, ofstream&); 
};

/****************************************************************************/

#endif


