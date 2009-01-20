/*
 *  readlistcommand.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 1/20/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#ifndef READLISTFILECOMMAND_H
#define READLISTFILECOMMAND_H
/*
 *  readlistcommand.h
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
#include "groupmap.h"
#include "sharedcommand.h"
#include "parselistcommand.h"


/* The read.list command parameter options are listfile and groupfile.  
The read.list command should be in the following format: 
read.shared(listfile=yourListFile, groupfile=yourGroupFile).  
The listfile parameter and groupfile paramaters are required. */		


class GlobalData;

class ReadListFileCommand : public Command {
public:
	ReadListFileCommand();
	~ReadListFileCommand();
	int execute();
	
private:
	GlobalData* globaldata;
	Command* shared;
	Command* parselist;
	GroupMap* groupMap;
	ReadMatrix* read;
	InputData* input;
	string filename;
};

#endif