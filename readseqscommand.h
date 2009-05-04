#ifndef READSEQSCOMMAND_H
#define READSEQSCOMMAND_H
/*
 *  readseqscommand.h
 *  Mothur
 *
 *  Created by Thomas Ryabin on 4/13/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "command.hpp"
#include "readfasta.h"
#include "readnexus.h"
#include "readclustal.h"
#include "readseqsphylip.h"
#include "inputdata.h"
#include "groupmap.h"
#include "sharedcommand.h"
#include "parselistcommand.h"

/* The read.otu must be run before you execute a collect.single, rarefaction.single, summary.single, 
collect.shared, rarefaction.shared or summary.shared command. Mothur will generate a .list, .rabund and .sabund 
upon completion of the cluster command or you may use your own. The read.otu command parameter options are 
listfile, rabundfile, sabundfile, groupfile and orderfile. The reaad.otu command can be used in two ways. 
The first is to read a listfile, rabundfile or sabundfile and run the collect.single, rarefaction.single or summary.single. 
For this use the read.otu command should be in the following format: read.otu(listfile=yourListFile, orderfile=yourOrderFile). 
The listfile, rabundfile or sabundfile parameter is required, but you may only use one of them. 
The second way to use the read.otu command is to read a listfile and a groupfile so you can use the collect.shared, 
rarefaction.shared or summary.shared commands. In this case the read.otu command should be in the following format: 
read.otu(listfile=yourListFile, groupfile=yourGroupFile). The listfile parameter and groupfile paramaters are required. 
When using the command the second way read.otu command parses the .list file and separates it into groups. 
It outputs a .shared file containing the OTU information for each group. The read.otu command also outputs a .list file for each group. */

class GlobalData;

class ReadSeqsCommand : public Command {
public:
	ReadSeqsCommand();
	~ReadSeqsCommand();
	int execute();
	
private:
	GlobalData* globaldata;
	ReadFasta* readFasta;
	ReadNexus* readNexus;
	ReadClustal* readClustal;
	ReadPhylip* readPhylip;
	InputData* input;
	Command* shared;
	Command* parselist;
	string filename;
};

#endif
