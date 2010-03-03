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
#include "sequence.hpp"

class TrimSeqsCommand : public Command {
public:
	TrimSeqsCommand(string);
	~TrimSeqsCommand();
	int execute();
	void help();
	
private:
	void getOligos(vector<ofstream*>&);
	bool stripQualThreshold(Sequence&, ifstream&);
	bool cullQualAverage(Sequence&, ifstream&);
	bool stripBarcode(Sequence&, int&);
	bool stripForward(Sequence&);
	bool stripReverse(Sequence&);
	bool cullLength(Sequence&);
	bool cullHomoP(Sequence&);
	bool cullAmbigs(Sequence&);
	bool compareDNASeq(string, string);

	bool abort;
	string fastaFile, oligoFile, qFileName, outputDir;
	
	bool flip, allFiles, qtrim;
	int numFPrimers, numRPrimers, maxAmbig, maxHomoP, minLength, maxLength, qThreshold, qAverage;
	vector<string> forPrimer, revPrimer, outputNames;
	map<string, int> barcodes;
	vector<string> groupVector;
};

#endif
