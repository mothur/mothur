#ifndef COLLECTCOMMAND_H
#define COLLECTCOMMAND_H
/*
 *  collectcommand.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "command.hpp"
#include "ordervector.hpp"
#include "inputdata.h"
#include "collect.h"
#include "display.h"
#include "validcalculator.h"

/*The collect() command:
	The collect command generates a collector's curve from the given file.  
	The collect command outputs a file for each estimator you choose to use.  The collect command parameters are label, freq, single, abund.  
	No parameters are required.  
	The collect command should be in the following format: collect(label=yourLabel, freq=yourFreq, single=yourEstimators, abund=yourAbund). 
	example collect(label=unique-.01-.03, freq=10, single=collect-chao-ace-jack).  
	The default values for  freq is 100, for abund is 10, and single are collect-chao-ace-jack-bootstrap-shannon-npshannon-simpson.  
	The valid single estimators are: collect-chao-ace-jack-bootstrap-shannon-npshannon-simpson. 
	The label parameter is used to analyze specific labels in your input. */


class CollectCommand : public Command {
	
public:
	CollectCommand(string);	
	CollectCommand();
	~CollectCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "collect.single";			}
	string getCommandCategory()		{ return "OTU-Based Approaches";	}
	string getCitation() { return "Schloss PD, Handelsman J (2006). Introducing SONS, A tool that compares the membership of microbial communities. Appl Environ Microbiol 72: 6773-9. \nhttp://www.mothur.org/wiki/Collect.single"; }
	string getHelpString();	
	string getDescription()		{ return "generates collector's curves using calculators, that describe the richness, diversity, and other features of individual samples"; }
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
	
private:
	OrderVector* order;
	InputData* input;
	Collect* cCurve;
	vector<Display*> cDisplays;
	int abund, size;
	float freq;
	vector<string> outputNames;

	bool abort, allLines;
	set<string> labels; //holds labels to be used
	string label, calc, outputDir, sharedfile, listfile, rabundfile, sabundfile, format, inputfile;
	vector<string>  Estimators;
	vector<string> inputFileNames;
	vector<string> groups;
	
	vector<string> parseSharedFile(string);


};

#endif
