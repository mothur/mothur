#ifndef SUMMARYCOMMAND_H
#define SUMMARYCOMMAND_H
/*
 *  summarycommand.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "command.hpp"
#include "sabundvector.hpp"
#include "inputdata.h"
#include "calculator.h"
#include "readotu.h"
#include "validcalculator.h"

class GlobalData;

class SummaryCommand : public Command {

public:
	SummaryCommand(string);
	~SummaryCommand();
	int execute();
	void help();
	
private:
	GlobalData* globaldata;
	ReadOTUFile* read;
	vector<Calculator*> sumCalculators;	
	InputData* input;
	ValidCalculators* validCalculator;
	SAbundVector* sabund;
	int abund, size;

	bool abort, allLines;
	set<string> labels; //holds labels to be used
	string label, calc, outputDir;
	vector<string>  Estimators;
	vector<string> inputFileNames;
	vector<string> groups;
	
	vector<string> parseSharedFile(string);


};
#endif
