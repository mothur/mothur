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
#include "qualityscores.h"

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

	void getOligos(vector<string>&, vector<string>&);
	int stripBarcode(Sequence&, QualityScores&, int&);
	int stripForward(Sequence&, QualityScores&, int&);
	bool stripReverse(Sequence&, QualityScores&);
	bool cullLength(Sequence&);
	bool cullHomoP(Sequence&);
	bool cullAmbigs(Sequence&);
	bool compareDNASeq(string, string);
	int countDiffs(string, string);

	bool abort;
	string fastaFile, oligoFile, qFileName, outputDir;
	
	bool flip, allFiles, qtrim;
	int numFPrimers, numRPrimers, maxAmbig, maxHomoP, minLength, maxLength, processors, tdiffs, bdiffs, pdiffs, comboStarts;
	int qWindowSize, qWindowStep;
	double qRollAverage, qThreshold, qWindowAverage, qAverage;
	vector<string> revPrimer, outputNames;
	set<string> filesToRemove;
	map<string, int> barcodes;
	vector<string> groupVector;
	map<string, int> primers;
	map<string, int> combos;
	
	vector<int> processIDS;   //processid
	vector<linePair*> lines;
	vector<linePair*> qLines;
	
	int driverCreateTrim(string, string, string, string, string, string, string, vector<string>, vector<string>, linePair*, linePair*);	
	int createProcessesCreateTrim(string, string, string, string, string, string, string, vector<string>, vector<string>);
	int setLines(string, vector<linePair*>&);
	
};

#endif
