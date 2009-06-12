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
	~GetRAbundCommand();
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
	OptionParser* parser;
	map<string, string> parameters;
	map<string, string>::iterator it;
	bool abort, allLines;
	set<int> lines; //hold lines to be used
	set<string> labels; //holds labels to be used
	string line, label;

	
};

#endif

