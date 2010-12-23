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
#include "groupmap.h"

class TrimSeqsCommand : public Command {
public:
	TrimSeqsCommand(string);
	TrimSeqsCommand();
	~TrimSeqsCommand();
	vector<string> getRequiredParameters();
	vector<string> getValidParameters();
	vector<string> getRequiredFiles();
	map<string, vector<string> > getOutputFiles() { return outputTypes; }
	int execute();
	void help();
	
private:
	
	GroupMap* groupMap;
	
	struct linePair {
		unsigned long int start;
		unsigned long int end;
		linePair(unsigned long int i, unsigned long int j) : start(i), end(j) {}
	};
	
	void getOligos(vector<string>&, vector<string>&);
	int stripBarcode(Sequence&, QualityScores&, int&);
	int stripForward(Sequence&, QualityScores&, int&);
	bool stripReverse(Sequence&, QualityScores&);
	
	bool keepFirstTrim(Sequence&, QualityScores&);
	bool removeLastTrim(Sequence&, QualityScores&);

	bool cullLength(Sequence&);
	bool cullHomoP(Sequence&);
	bool cullAmbigs(Sequence&);
	bool compareDNASeq(string, string);
	int countDiffs(string, string);
	map<string, vector<string> > outputTypes;

	bool abort;
	string fastaFile, oligoFile, qFileName, groupfile, outputDir;
	
	bool flip, allFiles, qtrim;
	int numFPrimers, numRPrimers, maxAmbig, maxHomoP, minLength, maxLength, processors, tdiffs, bdiffs, pdiffs, comboStarts;
	int qWindowSize, qWindowStep, keepFirst, removeLast;
	double qRollAverage, qThreshold, qWindowAverage, qAverage;
	vector<string> revPrimer, outputNames;
	set<string> filesToRemove;
	map<string, int> barcodes;
	vector<string> groupVector;
	map<string, int> primers;
	map<string, int> combos;
	map<string, int> groupToIndex;
	
	vector<int> processIDS;   //processid
	vector<linePair*> lines;
	vector<linePair*> qLines;
	
	int driverCreateTrim(string, string, string, string, string, string, string, vector<string>, vector<string>, linePair*, linePair*);	
	int createProcessesCreateTrim(string, string, string, string, string, string, string, vector<string>, vector<string>);
	int setLines(string, string, vector<unsigned long int>&, vector<unsigned long int>&);
};

#endif
