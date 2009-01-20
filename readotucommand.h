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

#include <Carbon/Carbon.h>
#include <iostream>
#include <fstream>
#include "command.hpp"
#include "readmatrix.hpp"
#include "inputdata.h"

/* The read.otu command parameter options are listfile, sabundfile, rabundfile and orderfile.  
The read.otu command should be in the following format: 
read.otu(listfile=yourListFile, orderfile=yourOrderFile). The listfile, sabundfile or rabundfile are required, but only one may be used. */

class GlobalData;

class ReadOtuCommand : public Command {
public:
	ReadOtuCommand();
	~ReadOtuCommand();
	int execute();
	
private:
	GlobalData* globaldata;
	ReadMatrix* read;
	InputData* input;
	string filename;
};

#endif