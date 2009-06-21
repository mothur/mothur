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
#include "readotu.h"
#include "sharedlistvector.h"
#include "venn.h"
#include "validcalculator.h"


class GlobalData;


class VennCommand : public Command {

public:
	VennCommand(string);
	~VennCommand();
	int execute();
	void help();
	
private:
	GlobalData* globaldata;
	ReadOTUFile* read;
	InputData* input;
	SharedListVector* SharedList;
	Venn* venn;
	vector<Calculator*> vennCalculators;	
	ValidCalculators* validCalculator;
	vector<SharedRAbundVector*> lookup;
	SAbundVector* sabund;
	int abund;
	
	bool abort, allLines;
	set<int> lines; //hold lines to be used
	set<string> labels; //holds labels to be used
	string format, groups, calc, line, label;
	vector<string> Estimators, Groups;


};



#endif
