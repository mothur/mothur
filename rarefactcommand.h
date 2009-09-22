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

#include "command.hpp"
#include "ordervector.hpp"
#include "inputdata.h"
#include "rarefact.h"
#include "display.h"
#include "readotu.h"
#include "validcalculator.h"


class GlobalData;

class RareFactCommand : public Command {
	
public:
	RareFactCommand(string);	
	~RareFactCommand();
	int execute();
	void help();	
	
private:
	GlobalData* globaldata;
	vector<Display*> rDisplays;
	ReadOTUFile* read;
	OrderVector* order;
	InputData* input;
	ValidCalculators* validCalculator;
	Rarefact* rCurve;
	int freq, nIters, abund;

	bool abort, allLines;
	set<string> labels; //holds labels to be used
	string label, calc;
	vector<string>  Estimators;

};

#endif
