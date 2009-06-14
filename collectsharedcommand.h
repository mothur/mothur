#ifndef COLLECTSHAREDCOMMAND_H
#define COLLECTSHAREDCOMMAND_H
/*
 *  collectsharedcommand.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */


#include "command.hpp"
#include "sharedordervector.h"
#include "sharedlistvector.h"
#include "inputdata.h"
#include "groupmap.h"
#include "collect.h"
#include "display.h"
#include "readotu.h"
#include "validcalculator.h"
#include "sharedutilities.h"

/* The collect.shared() command:
	The collect command generates a collector's curve from the given file representing several groups.  
	The collect.shared command can only be executed after a successful read.shared command. 
	It outputs a file for each estimator you choose to use.  The collect.shared command parameters are label, line, freq and shared.  
	No parameters are required, but you may not use both the line and label parameters at the same time. 
	The collect.shared command should be in the following format: collect.shared(label=yourLabel, line=yourLines, 
	freq=yourFreq, shared=yourEstimators). Example collect.shared(label=unique-.01-.03, line=0,5,10, freq=10, 
	shared=sharedChao-sharedAce-sharedJabund). The default value for
	freq is 100 and shared are sharedChao-sharedAce-sharedJabund-sharedSorensonAbund-sharedJclass-sharedSorClass-sharedJest-sharedSorEst-SharedThetaYC-SharedThetaN.  
	The valid shared estimators are: sharedChao-sharedAce-sharedJabund-sharedSorensonAbund-sharedJclass-sharedSorClass-sharedJest-sharedSorEst-SharedThetaYC-SharedThetaN.  
	The label and line parameters are used to analyze specific lines in your input. */


class GlobalData;

class CollectSharedCommand : public Command {
	
public:
	CollectSharedCommand(string);	
	~CollectSharedCommand();
	int execute();	
	void help();
	
private:
	GlobalData* globaldata;
	SharedUtil* util;
	GroupMap* groupmap;
	SharedListVector* SharedList;
	ReadOTUFile* read;
	SharedOrderVector* order;
	InputData* input;
	ValidCalculators* validCalculator;
	Collect* cCurve;
	vector<Display*> cDisplays;
	int freq;
	string format;

	bool abort, allLines;
	set<int> lines; //hold lines to be used
	set<string> labels; //holds labels to be used
	string line, label, calc, groups;
	vector<string>  Estimators, Groups;


};

#endif
