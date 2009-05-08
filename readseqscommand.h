#ifndef READSEQSCOMMAND_H
#define READSEQSCOMMAND_H
/*
 *  readseqscommand.h
 *  Mothur
 *
 *  Created by Thomas Ryabin on 4/13/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "command.hpp"
#include "readfasta.h"
#include "readnexus.h"
#include "readclustal.h"
#include "readseqsphylip.h"

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
	string filename;
};

#endif
