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
#include "readotu.h"
#include "listvector.hpp"

class GlobalData;

class GetRAbundCommand : public Command {
public:
	GetRAbundCommand(string);
	GetRAbundCommand();
	~GetRAbundCommand();
	vector<string> getRequiredParameters();
	vector<string> getValidParameters();
	vector<string> getRequiredFiles();
	map<string, vector<string> > getOutputFiles() { return outputTypes; }
	int execute();
	void help();
	
private:
	GlobalData* globaldata;
	string filename;
	ofstream out;
	ReadOTUFile* read;
	InputData* input;
	ListVector* list;
	RAbundVector* rabund;
	vector<string> outputNames;
	map<string, vector<string> > outputTypes;

	bool abort, allLines, sorted;
	set<string> labels; //holds labels to be used
	string label;

	
};

#endif

