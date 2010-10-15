#ifndef SUMMARYSHAREDCOMMAND_H
#define SUMMARYSHAREDCOMMAND_H
/*
 *  summarysharedcommand.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */


#include "command.hpp"
#include "sharedrabundvector.h"
#include "inputdata.h"
#include "calculator.h"
#include "readotu.h"
#include "validcalculator.h"

class GlobalData;

class SummarySharedCommand : public Command {

public:
	SummarySharedCommand(string);
	SummarySharedCommand();
	~SummarySharedCommand();
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
	vector<linePair> lines;
	GlobalData* globaldata;
	ReadOTUFile* read;
	vector<Calculator*> sumCalculators;	
	InputData* input;
	ValidCalculators* validCalculator;
	
	bool abort, allLines, mult, all, createPhylip;
	set<string> labels; //holds labels to be used
	string label, calc, groups;
	vector<string>  Estimators, Groups, outputNames;
	map<string, vector<string> > outputTypes;
	vector<SharedRAbundVector*> lookup;
	string format, outputDir;
	int numGroups, processors;
	int process(vector<SharedRAbundVector*>, string, string);
	int driver(vector<SharedRAbundVector*>, int, int, string, string, vector< vector<seqDist> >&);

};

#endif
