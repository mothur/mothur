#ifndef HEATMAPSIMCOMMAND_H
#define HEATMAPSIMCOMMAND_H

/*
 *  heatmapsimcommand.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 6/8/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "command.hpp"
#include "inputdata.h"
#include "readotu.h"
#include "validcalculator.h"
#include "heatmapsim.h"

class GlobalData;

class HeatMapSimCommand : public Command {

public:
	HeatMapSimCommand(string);
	~HeatMapSimCommand();
	int execute();
	void help();
	
private:
	GlobalData* globaldata;
	ReadOTUFile* read;
	InputData* input;
	vector<SharedRAbundVector*> lookup;
	vector<Calculator*> heatCalculators;
	ValidCalculators* validCalculator;
	HeatMapSim* heatmap;
	OptionParser* parser;
	map<string, string> parameters;
	map<string, string>::iterator it;
	bool abort, allLines;
	set<int> lines; //hold lines to be used
	set<string> labels; //holds labels to be used
	string format, groups, line, label, calc;
	vector<string> Estimators, Groups;


};

#endif

