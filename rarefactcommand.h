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
#include "validcalculator.h"

class RareFactCommand : public Command {
	
public:
	RareFactCommand(string);
	RareFactCommand();	
	~RareFactCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "rarefaction.single";		}
	string getCommandCategory()		{ return "OTU-Based Approaches";	}
	string getHelpString();	
	string getCitation() { return "Magurran AE (2004). Measuring biological diversity. Blackwell Pub.: Malden, Ma. \nhttp://www.mothur.org/wiki/Rarefaction.single"; }
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	
	vector<Display*> rDisplays;
	OrderVector* order;
	InputData* input;
	Rarefact* rCurve;
	int nIters, abund, processors;
	float freq;
	
	bool abort, allLines;
	set<string> labels; //holds labels to be used
	string label, calc, sharedfile, listfile, rabundfile, sabundfile, format, inputfile;
	vector<string>  Estimators;
	vector<string> inputFileNames, outputNames;
	vector<string> groups;
	string outputDir;
	
	vector<string> parseSharedFile(string);
};

#endif
