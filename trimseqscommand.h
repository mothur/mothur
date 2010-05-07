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

	struct linePair {
		int start;
		int num;
		linePair(long int i, int j) : start(i), num(j) {}
	};

	void getOligos(vector<string>&);
	bool stripQualThreshold(Sequence&, ifstream&);
	bool cullQualAverage(Sequence&, ifstream&);
	bool stripBarcode(Sequence&, int&);
	bool stripForward(Sequence&);
	bool stripReverse(Sequence&);
	bool cullLength(Sequence&);
	bool cullHomoP(Sequence&);
	bool cullAmbigs(Sequence&);
	bool compareDNASeq(string, string);
	bool compareDNASeq(string, string, int, int&);

	bool abort;
	string fastaFile, oligoFile, qFileName, outputDir;
	
	bool flip, allFiles, qtrim;
	int numFPrimers, numRPrimers, maxAmbig, maxHomoP, minLength, maxLength, qThreshold, qAverage, processors, diffs;
	vector<string> forPrimer, revPrimer, outputNames;
	map<string, int> barcodes;
	vector<string> groupVector;
	
	vector<int> processIDS;   //processid
	vector<linePair*> lines;
	vector<linePair*> qLines;
	
	int driverCreateTrim(string, string, string, string, string, vector<string>, linePair*, linePair*);	
	int createProcessesCreateTrim(string, string, string, string, string, vector<string>);
	int setLines(string, vector<linePair*>&);
	
};

#endif
