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
	~GetSAbundCommand();
	int execute();
	void help();
	
private:
	GlobalData* globaldata;
	string filename, format;
	ofstream out;
	ReadOTUFile* read;
	OrderVector* order;
	OrderVector* lastOrder;
	InputData* input;
	SAbundVector* sabund;

	bool abort, allLines;
	set<int> lines; //hold lines to be used
	set<string> labels; //holds labels to be used
	string line, label;

};

#endif
