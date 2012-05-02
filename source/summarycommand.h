#ifndef SUMMARYCOMMAND_H
#define SUMMARYCOMMAND_H
/*
 *  summarycommand.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "command.hpp"
#include "sabundvector.hpp"
#include "inputdata.h"
#include "calculator.h"
#include "validcalculator.h"

class SummaryCommand : public Command {

public:
	SummaryCommand(string);
	SummaryCommand();
	~SummaryCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "summary.single";			}
	string getCommandCategory()		{ return "OTU-Based Approaches";	}
	string getHelpString();	
	string getCitation() { return "http://www.mothur.org/wiki/Summary.single"; }
	string getDescription()		{ return "generate summary file that has the calculator value for each line in the OTU data"; }

	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	vector<Calculator*> sumCalculators;	
	InputData* input;
	SAbundVector* sabund;
	int abund, size;

	bool abort, allLines, groupMode;
	set<string> labels; //holds labels to be used
	string label, calc, outputDir, sharedfile, listfile, rabundfile, sabundfile, format, inputfile;
	vector<string>  Estimators;
	vector<string> inputFileNames, outputNames;
	vector<string> groups;
	
	vector<string> parseSharedFile(string);
	string createGroupSummaryFile(int, int, vector<string>&);


};
#endif
