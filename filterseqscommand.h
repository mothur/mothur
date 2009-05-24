#ifndef FILTERSEQSCOMMAND_H
#define FILTERSEQSCOMMAND_H

/*
 *  filterseqscommand.h
 *  Mothur
 *
 *  Created by Thomas Ryabin on 5/4/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "command.hpp"
#include "mothur.h"
#include "globaldata.hpp"
#include "readfasta.h"
#include "readnexus.h"
#include "readclustal.h"
#include "readseqsphylip.h"

using namespace std;

class FilterSeqsCommand : public Command {

public:
	FilterSeqsCommand();
	~FilterSeqsCommand() {};
	int execute();	
	
private:
	void doTrump();
	void doSoft();
	void doHard();
	void doVertical();
	
	int alignmentLength;
	int numSeqs;
	
	GlobalData* globaldata;	
	ReadSeqs* readSeqs;
	SequenceDB* db;
	
	string filter;

};

#endif