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
	~SummarySharedCommand();
	int execute();
	void help();
	
private:
	GlobalData* globaldata;
	ReadOTUFile* read;
	vector<Calculator*> sumCalculators;	
	InputData* input;
	ValidCalculators* validCalculator;
	
	bool abort, allLines, mult, all;
	set<string> labels; //holds labels to be used
	string label, calc, groups;
	vector<string>  Estimators, Groups;
	vector<SharedRAbundVector*> lookup;
	string outputFileName, format, outAllFileName, outputDir;
	ofstream outputFileHandle, outAll;
	void process(vector<SharedRAbundVector*>);

};

#endif
