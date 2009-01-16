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

#include <Carbon/Carbon.h>
#include <iostream>
#include <fstream>
#include <vector>
#include "command.hpp"
#include "sharedordervector.h"
#include "listvector.hpp"
#include "inputdata.h"
#include "groupmap.h"
#include "rarefact.h"
#include "display.h"
#include "readmatrix.hpp"

/* The rarefaction.shared() command:
	The rarefaction command generates a rarefaction curve from a given file representing several groups.  
	The rarefaction.shared command can only be executed after a successful read.shared command. It outputs a file for each estimator you choose to use.  
	The rarefaction.shared command parameters are label, line, iters, jumble and sharedrarefaction.  
	No parameters are required, but you may not use both the line and label  parameters at the same time.  
	The rarefaction command should be in the following format: rarefaction.shared(label=yourLabel, line=yourLines, iters=yourIters, 
	jumble= yourJumble, sharedrarefaction=yourEstimators).  Example rarefaction.shared(label=unique-.01-.03, line=0,5,10, iters=10000, 
	jumble=1, sharedrarefaction =sharedobserved).  The default values for jumble is 0 (meaning don’t jumble, if it’s set to 1 then it will jumble), 
	iters is 1000 and sharedrarefaction is sharedobserved which calculates the shared rarefaction curve for the observed richness. 
	 The valid sharedrarefaction estimator is sharedobserved. The label and line parameters are used to analyze specific lines in your input. */


class GlobalData;

class RareFactSharedCommand : public Command {
	
public:
	RareFactSharedCommand();	
	~RareFactSharedCommand();
	int execute();	
	
private:
	GlobalData* globaldata;
	GroupMap* groupmap;
	ListVector* list;
	ReadMatrix* read;
	SharedOrderVector* order;
	InputData* input;
	Rarefact* rCurve;
	vector<Display*> rDisplays;
	int freq, nIters;

};

#endif