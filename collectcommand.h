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
#include "groupmap.h"
#include "collect.h"
#include "display.h"
#include "readotu.h"
#include "validcalculator.h"

/*The collect() command:
	The collect command generates a collector's curve from the given file.  
	The collect command can only be executed after a successful read.list, read.sabund or read.rabund command, with one exception. 
	The collect command can be executed after a successful cluster command.  It will use the .list file from the output of the cluster.  
	The collect command outputs a file for each estimator you choose to use.  The collect command parameters are label, freq, single, abund.  
	No parameters are required.  
	The collect command should be in the following format: collect(label=yourLabel, freq=yourFreq, single=yourEstimators, abund=yourAbund). 
	example collect(label=unique-.01-.03, freq=10, single=collect-chao-ace-jack).  
	The default values for  freq is 100, for abund is 10, and single are collect-chao-ace-jack-bootstrap-shannon-npshannon-simpson.  
	The valid single estimators are: collect-chao-ace-jack-bootstrap-shannon-npshannon-simpson. 
	The label parameter is used to analyze specific labels in your input. */



class GlobalData;

class CollectCommand : public Command {
	
public:
	CollectCommand(string);	
	~CollectCommand();
	int execute();	
	void help();
	
private:
	GlobalData* globaldata;
	ReadOTUFile* read;
	OrderVector* order;
	InputData* input;
	Collect* cCurve;
	ValidCalculators* validCalculator;
	vector<Display*> cDisplays;
	int freq, abund, size;

	bool abort, allLines;
	set<string> labels; //holds labels to be used
	string label, calc, outputDir;
	vector<string>  Estimators;
	vector<string> inputFileNames;
	vector<string> groups;
	
	vector<string> parseSharedFile(string);


};

#endif
