#ifndef HEATMAPCOMMAND_H
#define HEATMAPCOMMAND_H

/*
 *  heatmapcommand.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 3/25/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */


#include "command.hpp"
#include "inputdata.h"
#include "readotu.h"
#include "sharedlistvector.h"
#include "heatmap.h"
#include "rabundvector.hpp"

class GlobalData;

class HeatMapCommand : public Command {

public:
	HeatMapCommand(string);
	~HeatMapCommand();
	int execute();
	void help();
	
private:
	GlobalData* globaldata;
	ReadOTUFile* read;
	InputData* input;
	SharedListVector* SharedList;
	RAbundVector* rabund;
	vector<SharedRAbundVector*> lookup;
	HeatMap* heatmap;

	bool abort, allLines;
	set<string> labels; //holds labels to be used
	string format, groups, sorted, scale, label, outputDir;
	vector<string> Groups;


};

#endif

