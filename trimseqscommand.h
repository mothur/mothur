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
	TrimSeqsCommand();
	~TrimSeqsCommand();
	int execute();
	
private:
	void getOligos();
	bool stripBarcode(Sequence&, string&);
	bool stripForward(Sequence&);
	bool stripReverse(Sequence&);
	bool cullLength(Sequence&);
	bool cullHomoP(Sequence&);
	bool cullAmbigs(Sequence&);
	
	GlobalData* globaldata;

	bool oligos, flip;
	int numFPrimers, numRPrimers, maxAmbig, maxHomoP, minLength, maxLength;
	vector<string> forPrimer, revPrimer;
	map<string, string> barcodes;
};

#endif
