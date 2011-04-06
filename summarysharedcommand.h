#ifndef SUMMARYSHAREDCOMMAND_H
#define SUMMARYSHAREDCOMMAND_H
/*
 *  summarysharedcommand.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */


#include "command.hpp"
#include "sharedrabundvector.h"
#include "inputdata.h"
#include "calculator.h"
#include "validcalculator.h"

class SummarySharedCommand : public Command {

public:
	SummarySharedCommand(string);
	SummarySharedCommand();
	~SummarySharedCommand() {}
	
	vector<string> setParameters();
	string getCommandName()			{ return "summary.shared";			}
	string getCommandCategory()		{ return "OTU-Based Approaches";	}
	string getHelpString();	
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
	
private:
	struct linePair {
		int start;
		int end;
	};
	vector<linePair> lines;
	vector<Calculator*> sumCalculators;	
	InputData* input;
	
	bool abort, allLines, mult, all, createPhylip;
	set<string> labels; //holds labels to be used
	string label, calc, groups, sharedfile;
	vector<string>  Estimators, Groups, outputNames;
	vector<SharedRAbundVector*> lookup;
	string format, outputDir;
	int numGroups, processors;
	int process(vector<SharedRAbundVector*>, string, string);
	int driver(vector<SharedRAbundVector*>, int, int, string, string, vector< vector<seqDist> >&);

};

#endif
