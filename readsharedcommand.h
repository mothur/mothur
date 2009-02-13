#ifndef READSHAREDCOMMAND_H
#define READSHAREDCOMMAND_H
/*
 *  readsharedcommand.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 2/12/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */
#include <iostream>
#include <fstream>
#include "command.hpp"
#include "readmatrix.hpp"
#include "inputdata.h"

/* The read.shared must be run before you execute a collect.shared, rarefaction.shared or summary.shared command.
The read.shared command is used to read a shared file. The read.shared should be entered in the following format:
read.shared(shared=yourSharedFile). The shared parameter is required.  */

class GlobalData;

class ReadSharedCommand : public Command {
public:
	ReadSharedCommand();
	~ReadSharedCommand();
	int execute();
	
private:
	GlobalData* globaldata;
	ReadMatrix* read;
	InputData* input;
	string filename;
};

#endif




