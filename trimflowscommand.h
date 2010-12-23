#ifndef TRIMFLOWSCOMMAND_H
#define TRIMFLOWSCOMMAND_H

/*
 *  trimflowscommand.h
 *  Mothur
 *
 *  Created by Pat Schloss on 12/22/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "mothur.h"
#include "command.hpp"
#include "sequence.hpp"
#include "flowdata.h"
#include "groupmap.h"

class TrimFlowsCommand : public Command {
public:
	TrimFlowsCommand(string);
	TrimFlowsCommand();
	~TrimFlowsCommand();
	vector<string> getRequiredParameters();
	vector<string> getValidParameters();
	vector<string> getRequiredFiles();
	map<string, vector<string> > getOutputFiles() { return outputTypes; }
	int execute();
	void help();
	
private:
	bool abort;

//	GroupMap* groupMap;
	
	struct linePair {
		unsigned long int start;
		unsigned long int end;
		linePair(unsigned long int i, unsigned long int j) : start(i), end(j) {}
	};
	int comboStarts;
	vector<int> processIDS;   //processid
	vector<linePair*> lines;
	vector<linePair*> qLines;
	map<string, vector<string> > outputTypes;
	vector<string> outputNames;
	set<string> filesToRemove;

	
	
	void getOligos(vector<vector<string> >&);		//a rewrite of what is in trimseqscommand.h
	int stripBarcode(Sequence&, int&);				//largely redundant with trimseqscommand.h
	int stripForward(Sequence&, int&);				//largely redundant with trimseqscommand.h
	bool stripReverse(Sequence&);					//largely redundant with trimseqscommand.h
	bool compareDNASeq(string, string);				//largely redundant with trimseqscommand.h
	int countDiffs(string, string);					//largely redundant with trimseqscommand.h

	
	bool allFiles;
//	int processors;
	int numFPrimers, numRPrimers;
	int totalFlows, minFlows, minLength, maxLength, maxHomoP, tdiffs, bdiffs, pdiffs;
	float signal, noise;
	bool fasta;
	
	
	string flowFileName, oligoFileName, outputDir;


	map<string, int> barcodes;
	map<string, int> primers;
	vector<string> revPrimer;

	vector<string> primerNameVector;	//needed here?
	vector<string> barcodeNameVector;	//needed here?

	map<string, int> combos;			//needed here?
	map<string, int> groupToIndex;		//needed here?
	
	
	int driverCreateTrim(string, string, string, string);
	
//	int createProcessesCreateTrim(string, string, string, string, string, string, string, vector<string>, vector<string>){};
	int setLines(string, string, vector<unsigned long int>&, vector<unsigned long int>&){};
	
};


#endif
