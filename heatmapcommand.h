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
#include "sharedlistvector.h"
#include "heatmap.h"
#include "rabundvector.hpp"


class HeatMapCommand : public Command {

public:
	HeatMapCommand(string);
	HeatMapCommand();
	~HeatMapCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "heatmap.bin";				}
	string getCommandCategory()		{ return "OTU-Based Approaches";	}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "http://www.mothur.org/wiki/Heatmap.bin"; }
	string getDescription()		{ return "generate a heatmap where the color represents the relative abundanceof an OTU"; }

	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
	
private:
	InputData* input;
	RAbundVector* rabund;
	vector<SharedRAbundVector*> lookup;
	vector<SharedRAbundFloatVector*> lookupFloat;
	HeatMap* heatmap;

	bool abort, allLines;
	set<string> labels; //holds labels to be used
	string format, groups, sorted, scale, label, outputDir, sharedfile, relabundfile, listfile, rabundfile, sabundfile, inputfile;
	vector<string> Groups, outputNames;
	int numOTU, fontSize;


};

#endif

