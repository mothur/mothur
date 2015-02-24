#ifndef VENNCOMMAND_H
#define VENNCOMMAND_H
/*
 *  venncommand.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 3/30/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */
 
#include "command.hpp"
#include "inputdata.h"
#include "sharedlistvector.h"
#include "venn.h"
#include "validcalculator.h"

class VennCommand : public Command {

public:
	VennCommand(string);
	VennCommand();
	~VennCommand() {}
	
	vector<string> setParameters();
	string getCommandName()			{ return "venn";					}
	string getCommandCategory()		{ return "OTU-Based Approaches";	}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "http://www.mothur.org/wiki/Venn"; }
	string getDescription()		{ return "generates a Venn diagram from data provided in a shared file"; }

	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	InputData* input;
	SharedListVector* SharedList;
	Venn* venn;
	vector<Calculator*> vennCalculators;	
	vector<SharedRAbundVector*> lookup;
	set< set<int> > combos;
	SAbundVector* sabund;
	int abund, fontsize, perm;
	
	bool abort, allLines, nseqs, sharedOtus;
	set<string> labels; //holds labels to be used
	string format, groups, calc, label, outputDir, sharedfile, listfile, inputfile;
	vector<string> Estimators, Groups, outputNames;
	
	set< set<int> > findCombinations(int);
	int getCombos(set<int>, set< set<int> >&);


};



#endif
