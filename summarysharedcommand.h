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


/*The summary.shared() command
	The summary.shared command can only be executed after a successful read.shared command. 
	It outputs a file for each estimator you choose to use.  The summary.shared command parameters are label, 
	line and sharedsummary.  No parameters are required, but you may not use both the line and label parameters at the same time.  
	The summary.shared command should be in the following format: summary.shared(label=yourLabel, 
	line=yourLines, sharedsummary=yourEstimators).  
	Example summary.shared(label=unique-.01-.03, line=0,5,10, sharedsummary=sharedChao-sharedAce-sharedJabund
	-sharedSorensonAbund-sharedJclass-sharedSorClass-sharedJest-sharedSorEst-SharedThetaYC-SharedThetaN).  
	The default value for sharedsummary is sharedChao-sharedAce-sharedJabund-sharedSorensonAbund-sharedJclass-sharedSorClass-sharedJest-sharedSorEst-SharedThetaYC-SharedThetaN. 
	The valid sharedsummary estimators are: sharedChao-sharedAce-sharedJabund-sharedSorensonAbund-sharedJclass-sharedSorClass
	-sharedJest-sharedSorEst-SharedThetaYC-SharedThetaN.  The label and line parameters are used to analyze specific lines in your input. */


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
	SharedListVector* SharedList;
	OptionParser* parser;
	map<string, string> parameters;
	map<string, string>::iterator it;
	bool abort, allLines, mult;
	set<int> lines; //hold lines to be used
	set<string> labels; //holds labels to be used
	string line, label, calc, groups;
	vector<string>  Estimators, Groups;
	vector<SharedRAbundVector*> lookup;
	string outputFileName, format, outAllFileName;
	ofstream outputFileHandle, outAll;
	void process(vector<SharedRAbundVector*>);

};

#endif
