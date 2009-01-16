#ifndef READLISTFILECOMMAND_H
#define READLISTFILECOMMAND_H
/*
 *  readlistfilecommand.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include <Carbon/Carbon.h>
#include <iostream>
#include <fstream>
#include "command.hpp"
#include "readmatrix.hpp"
#include "inputdata.h"

/* The read.list command parameter options are listfile and orderfile.  
The read.list command should be in the following format: 
read.list(listfile=yourListFile, orderfile=yourOrderFile). The listfile parameter is required. */

class GlobalData;

class ReadListFileCommand : public Command {
public:
	ReadListFileCommand();
	~ReadListFileCommand();
	int execute();
	
private:
	GlobalData* globaldata;
	ReadMatrix* read;
	InputData* input;
	string filename;
};

#endif