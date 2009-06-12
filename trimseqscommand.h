#ifndef TRIMSEQSCOMMAND_H
#define TRIMSEQSCOMMAND_H

/*
 *  trimseqscommand.h
 *  Mothur
 *
 *  Created by Pat Schloss on 6/6/09.
 *  Copyright 2009 Patrick D. Schloss. All rights reserved.
 *
 */

#include "mothur.h"
#include "command.hpp"
#include "globaldata.hpp"
#include "sequence.hpp"

class TrimSeqsCommand : public Command {
public:
	TrimSeqsCommand(string);
	~TrimSeqsCommand();
	int execute();
	void help();
	
private:
	void getOligos();
	bool stripBarcode(Sequence&, string&);
	bool stripForward(Sequence&);
	bool stripReverse(Sequence&);
	bool cullLength(Sequence&);
	bool cullHomoP(Sequence&);
	bool cullAmbigs(Sequence&);
	
	GlobalData* globaldata;
	OptionParser* parser;
	map<string, string> parameters;
	map<string, string>::iterator it;
	bool abort;
	string fastafile;
	bool oligos, flip;
	int numFPrimers, numRPrimers, maxAmbig, maxHomoP, minLength, maxLength;
	vector<string> forPrimer, revPrimer;
	map<string, string> barcodes;
};

#endif
