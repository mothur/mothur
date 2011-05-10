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
#include "inputdata.h"
#include "collect.h"
#include "display.h"
#include "validcalculator.h"
#include "sharedutilities.h"


class CollectSharedCommand : public Command {
	
public:
	CollectSharedCommand(string);	
	CollectSharedCommand();	
	~CollectSharedCommand();
	
	vector<string> setParameters();
	string getCommandName()			{ return "collect.shared";			}
	string getCommandCategory()		{ return "OTU-Based Approaches";	}
	string getHelpString();	
	string getCitation() { return "http://www.mothur.org/wiki/Collect.shared"; }
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
	
private:
	SharedOrderVector* order;
	InputData* input;
	Collect* cCurve;
	vector<Display*> cDisplays;
	float freq;

	bool abort, allLines, all;
	set<string> labels; //holds labels to be used
	string label, calc, groups, outputDir, sharedfile;
	vector<string>  Estimators, Groups, outputNames;
};

#endif
