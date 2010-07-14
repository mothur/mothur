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
		unsigned long int start;
		int num;
		linePair(unsigned long int i, int j) : start(i), num(j) {}
	};

	void getOligos(vector<string>&);
	bool stripQualThreshold(Sequence&, ifstream&);
	bool cullQualAverage(Sequence&, ifstream&);
	int stripBarcode(Sequence&, int&);
	int stripForward(Sequence&, int&);
	bool stripReverse(Sequence&);
	bool cullLength(Sequence&);
	bool cullHomoP(Sequence&);
	bool cullAmbigs(Sequence&);
	bool compareDNASeq(string, string);
	int countDiffs(string, string);//, int, int&, int);

	bool abort;
	string fastaFile, oligoFile, qFileName, outputDir;
	
	bool flip, allFiles, qtrim;
	int numFPrimers, numRPrimers, maxAmbig, maxHomoP, minLength, maxLength, qThreshold, qAverage, processors, tdiffs, bdiffs, pdiffs, comboStarts;
	vector<string> revPrimer, outputNames;
	set<string> filesToRemove;
	map<string, int> barcodes;
	vector<string> groupVector;
	map<string, int> primers;
	map<string, int> combos;
	
	vector<int> processIDS;   //processid
	vector<linePair*> lines;
	vector<linePair*> qLines;
	
	int driverCreateTrim(string, string, string, string, string, vector<string>, linePair*, linePair*);	
	int createProcessesCreateTrim(string, string, string, string, string, vector<string>);
	int setLines(string, vector<linePair*>&);
	
};

#endif
