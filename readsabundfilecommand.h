#ifndef READSABUNDFILECOMMAND_H
#define READSABUNDFILECOMMAND_H
/*
 *  readsabundfilecommand.h
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

/* The read.sabund command parameter options are sabundfile and orderfile.  
The read.sabund command should be in the following format: 
read.sabund(sabundfile=yourSabundFile, orderfile=yourOrderFile). The sabundfile parameter is required.*/

class GlobalData;

class ReadSAbundFileCommand : public Command {
public:
	ReadSAbundFileCommand();
	~ReadSAbundFileCommand();
	int execute();
	
private:
	GlobalData* globaldata;
	ReadMatrix* read;
	InputData* input;
	string filename;
};

#endif