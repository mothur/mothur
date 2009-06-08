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
	HeatMapSimCommand();
	~HeatMapSimCommand();
	int execute();
	
private:
	GlobalData* globaldata;
	ReadOTUFile* read;
	InputData* input;
	vector<SharedRAbundVector*> lookup;
	vector<Calculator*> heatCalculators;
	ValidCalculators* validCalculator;
	HeatMapSim* heatmap;

};

#endif

