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



class CollectSharedCommand : public Command {
	
public:
	CollectSharedCommand(string);	
	~CollectSharedCommand();
	
	vector<string> setParameters();
	string getCommandName()			{ return "collect.shared";			}
	string getCommandCategory()		{ return "OTU-Based Approaches";	}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "Schloss PD, Handelsman J (2006). Introducing SONS, A tool that compares the membership of microbial communities. Appl Environ Microbiol 72: 6773-9. \nhttp://www.mothur.org/wiki/Collect.shared"; }
	string getDescription()		{ return "generates collector's curves for calculators, which describe the similarity between communities or their shared richness"; }

	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
	
private:
	
	vector<Display*> cDisplays;
	float freq;

	bool abort, allLines, all;
	set<string> labels; //holds labels to be used
	string label, calc, groups, sharedfile;
	vector<string>  Estimators, Groups, outputNames;
};

#endif
