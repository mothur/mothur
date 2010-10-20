#ifndef PAIRWISESEQSCOMMAND_H
#define PAIRWISESEQSCOMMAND_H

/*
 *  pairwiseseqscommand.h
 *  Mothur
 *
 *  Created by westcott on 10/20/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "mothur.h"
#include "command.hpp"
#include "database.hpp"
#include "alignment.hpp"
#include "validcalculator.h"
#include "dist.h"
#include "sequencedb.h"

class PairwiseSeqsCommand : public Command {
	
public:
	PairwiseSeqsCommand(string);	
	PairwiseSeqsCommand();
	~PairwiseSeqsCommand();
	vector<string> getRequiredParameters();
	vector<string> getValidParameters();
	vector<string> getRequiredFiles();
	map<string, vector<string> > getOutputFiles() { return outputTypes; }
	int execute(); 
	void help();	
	
private:
	struct linePair {
		int start;
		int end;
	};
	map<int, int> processIDS;   //end line, processid
	vector<linePair*> lines;
	map<string, vector<string> > outputTypes;
	
	Alignment* alignment;
	Dist* distCalculator;
	SequenceDB alignDB;
	
	void createProcesses(string);
	int driver(int, int, string, float);
	int driver(int, int, string, string);
	
	#ifdef USE_MPI 
	int driverMPI(int, int, MPI_File&, float);
	int driverMPI(int, int, string, unsigned long int&);
	int driverMPI(int, int, string, unsigned long int&, string);
	#endif
	
	string fastaFileName, align, calc, outputDir, output;
	float match, misMatch, gapOpen, gapExtend, cutoff;
	int processors;
	vector<string> fastaFileNames, Estimators;
	vector<string> outputNames;
	
	bool abort, countends, compress;
};

#endif

