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
#include "globaldata.hpp"
#include "readfasta.h"
#include "readnexus.h"
#include "readclustal.h"
#include "readseqsphylip.h"


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
	string filter;	
	int alignmentLength;

	char trump;
	bool vertical;
	
	GlobalData* globaldata;	
//	ReadSeqs* readSeqs;
//	SequenceDB* db;
	

};

#endif
