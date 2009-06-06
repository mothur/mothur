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
	
	GlobalData* globaldata;
	
	int totalBarcodeCount, matchBarcodeCount;
	int totalFPrimerCount, matchFPrimerCount;
	int totalRPrimerCount, matchRPrimerCount;
	
	
	bool oligos, flip;
	int forwardPrimerMismatch, reversePrimerMismatch, barcodeMismatch;
	int numFPrimers, numRPrimers;
	vector<string> forPrimer, revPrimer;
	map<string, string> barcodes;

};

#endif
