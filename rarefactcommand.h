#ifndef RAREFACTCOMMAND_H
#define RAREFACTCOMMAND_H
/*
 *  rarefactcommand.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include <iostream>
#include <fstream>
#include <vector>
#include "command.hpp"
#include "ordervector.hpp"
#include "inputdata.h"
#include "rarefact.h"
#include "display.h"
#include "readmatrix.hpp"
#include "validcalculator.h"


/*The rarefaction() command:
	The rarefaction command generates a rarefaction curve from a given file.  
	The rarefaction command can only be executed after a successful read.list, read.sabund or read.rabund command, with one exception. 
	The rarefaction command can be executed after a successful cluster command.  It will use the .list file from the output of the cluster. 
	The rarefaction command outputs a file for each estimator you choose to use.  It is recommended to only use rarefaction estimator.  
	The rarefaction command parameters are label, line, iters, freq, rarefaction.  No parameters are required, 
	but you may not use both the line and label  parameters at the same time. The rarefaction command should be in the following format: 
	rarefaction(label=yourLabel, line=yourLines, iters=yourIters, freq=yourFreq, rarefaction=yourEstimators). 
	Example rarefaction(label=unique-.01-.03, line=0,5,10, iters=10000, freq=10, rarefaction=rarefaction-rchao-race-rjack-rbootstrap-rshannon-rnpshannon-rsimpson). 
	The default values for iters is 1000, freq is 100, and rarefaction is rarefaction which calculates the rarefaction curve for the observed richness. 
	The valid rarefaction estimators are: rarefaction-rchao-race-rjack-rbootstrap-rshannon-rnpshannon-rsimpson.  
	Rarefaction is the only recommended estimator.  The label and line parameters are used to analyze specific lines in your input. */
	

class GlobalData;

class RareFactCommand : public Command {
	
public:
	RareFactCommand();	
	~RareFactCommand();
	int execute();	
	
private:
	GlobalData* globaldata;
	vector<Display*> rDisplays;
	ReadMatrix* read;
	OrderVector* order;
	InputData* input;
	ValidCalculators* validCalculator;
	Rarefact* rCurve;
	int freq, nIters, abund;
};

#endif
