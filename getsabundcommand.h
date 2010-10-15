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
#include "readotu.h"
#include "sabundvector.hpp"

class GlobalData;

class GetSAbundCommand : public Command {
public:
	GetSAbundCommand(string);
	GetSAbundCommand();
	~GetSAbundCommand();
	vector<string> getRequiredParameters();
	vector<string> getValidParameters();
	vector<string> getRequiredFiles();
	map<string, vector<string> > getOutputFiles() { return outputTypes; }
	int execute();
	void help();
	
private:
	GlobalData* globaldata;
	string filename, format;
	ofstream out;
	ReadOTUFile* read;
	OrderVector* order;
	InputData* input;
	SAbundVector* sabund;
	vector<string> outputNames;
	map<string, vector<string> > outputTypes;

	bool abort, allLines;
	set<string> labels; //holds labels to be used
	string label;

};

#endif
