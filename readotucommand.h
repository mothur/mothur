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

class GlobalData;

class ReadOtuCommand : public Command {
public:
	ReadOtuCommand(string);
	~ReadOtuCommand();
	int execute();
	void help();
	
private:
	GlobalData* globaldata;
	InputData* input;
	Command* shared;
	GroupMap* groupMap;
	string filename, listfile, orderfile, sharedfile, label, groupfile, sabundfile, rabundfile, format;

	bool abort, allLines;
	set<string> labels; //holds labels to be used

};

#endif
