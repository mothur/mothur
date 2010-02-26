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

	bool abort, allLines, all;
	set<string> labels; //holds labels to be used
	string label, calc, groups, outputDir;
	vector<string>  Estimators, Groups, outputNames;


};

#endif
