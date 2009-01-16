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

#include <Carbon/Carbon.h>
#include <iostream>
#include <fstream>
#include <vector>
#include "command.hpp"
#include "sabundvector.hpp"
#include "inputdata.h"
#include "calculator.h"
#include "readmatrix.hpp"

/* The summary() command:
	The summary command can only be executed after a successful read.list, read.sabund or read.rabund command, with one exception. 
	The summary command can be executed after a successful cluster command.  It will use the .list file from the output of the cluster. 
	The summary command outputs a file for each estimator you choose to use.  The summary command parameters are label, line, summary. 
	 No parameters are required, but you may not use both the line and label  parameters at the same time. 
	 The summary command should be in the following format: summary(label=yourLabel, line=yourLines, summary=yourEstimators). 
	 Example summary(label=unique-.01-.03, line=0,5,10, summary=collect-chao-ace-jack-bootstrap-shannon-npshannon-simpson). 
	 The default value for summary is collect-chao-ace-jack-bootstrap-shannon-npshannon-simpson.  
	 The valid summary estimators are: collect-chao-ace-jack-bootstrap-shannon-npshannon-simpson.  
	 The label and line parameters are used to analyze specific lines in your input.  */


class GlobalData;

class SummaryCommand : public Command {

public:
	SummaryCommand();
	~SummaryCommand();
	int execute();
	
private:
	GlobalData* globaldata;
	ReadMatrix* read;
	vector<Calculator*> sumCalculators;	
	InputData* input;
	SAbundVector* sabund;
	string outputFileName;
	ofstream outputFileHandle;
};
#endif