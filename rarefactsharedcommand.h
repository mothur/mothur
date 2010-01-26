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

	bool abort, allLines, jumble;
	set<string> labels; //holds labels to be used
	string label, calc, groups, outputDir;
	vector<string>  Estimators, Groups;


};

#endif
