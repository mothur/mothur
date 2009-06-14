#ifndef RAREFACTSHAREDCOMMAND_H
#define RAREFACTSHAREDCOMMAND_H
/*
 *  rarefactsharedcommand.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/6/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "command.hpp"
#include "inputdata.h"
#include "rarefact.h"
#include "display.h"
#include "readotu.h"
#include "validcalculator.h"

/* The rarefaction.shared() command:
	The rarefaction command generates a rarefaction curve from a given file representing several groups.  
	The rarefaction.shared command can only be executed after a successful read.shared command. It outputs a file for each estimator you choose to use.  
	The rarefaction.shared command parameters are label, line, iters and sharedrarefaction.  
	No parameters are required, but you may not use both the line and label  parameters at the same time.  
	The rarefaction command should be in the following format: rarefaction.shared(label=yourLabel, line=yourLines, iters=yourIters, 
	 sharedrarefaction=yourEstimators).  Example rarefaction.shared(label=unique-.01-.03, line=0,5,10, iters=10000, 
	 sharedrarefaction =sharedobserved).  The default values for 
	iters is 1000 and sharedrarefaction is sharedobserved which calculates the shared rarefaction curve for the observed richness. 
	 The valid sharedrarefaction estimator is sharedobserved. The label and line parameters are used to analyze specific lines in your input. */


class GlobalData;

class RareFactSharedCommand : public Command {
	
public:
	RareFactSharedCommand(string);	
	~RareFactSharedCommand();
	int execute();	
	void help();
	
private:
	GlobalData* globaldata;
	ReadOTUFile* read;
	vector<SharedRAbundVector*> lookup;
	InputData* input;
	ValidCalculators* validCalculator;
	Rarefact* rCurve;
	vector<Display*> rDisplays;
	int freq, nIters;
	string format;

	bool abort, allLines;
	set<int> lines; //hold lines to be used
	set<string> labels; //holds labels to be used
	string line, label, calc, groups;
	vector<string>  Estimators, Groups;


};

#endif
