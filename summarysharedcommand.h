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
#include "readmatrix.hpp"
#include "sharedlistvector.h"
#include "validcalculator.h"
#include "sharedutilities.h"

/*The summary.shared() command
	The summary.shared command can only be executed after a successful read.shared command. 
	It outputs a file for each estimator you choose to use.  The summary.shared command parameters are label, 
	line, jumble and sharedsummary.  No parameters are required, but you may not use both the line and label parameters at the same time.  
	The summary.shared command should be in the following format: summary.shared(label=yourLabel, 
	line=yourLines, jumble=yourJumble, sharedsummary=yourEstimators).  
	Example summary.shared(label=unique-.01-.03, line=0,5,10, jumble=1, sharedsummary=sharedChao-sharedAce-sharedJabund
	-sharedSorensonAbund-sharedJclass-sharedSorClass-sharedJest-sharedSorEst-SharedThetaYC-SharedThetaN).  
	The default value for jumble is 0 (meaning don’t jumble, if it’s set to 1 then it will jumble) and 
	sharedsummary is sharedChao-sharedAce-sharedJabund-sharedSorensonAbund-sharedJclass-sharedSorClass-sharedJest-sharedSorEst-SharedThetaYC-SharedThetaN. 
	The valid sharedsummary estimators are: sharedChao-sharedAce-sharedJabund-sharedSorensonAbund-sharedJclass-sharedSorClass
	-sharedJest-sharedSorEst-SharedThetaYC-SharedThetaN.  The label and line parameters are used to analyze specific lines in your input. */


class GlobalData;


class SummarySharedCommand : public Command {

public:
	SummarySharedCommand();
	~SummarySharedCommand();
	int execute();
	
private:
	GlobalData* globaldata;
	ReadMatrix* read;
	SharedUtil* util;
	vector<Calculator*> sumCalculators;	
	InputData* input;
	ValidCalculators* validCalculator;
	SharedListVector* SharedList;
	SharedOrderVector* order;
	vector<SharedRAbundVector*> lookup;
	string outputFileName, format;
	ofstream outputFileHandle;

};

#endif
