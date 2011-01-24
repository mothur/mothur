#ifndef NORMALIZESHAREDCOMMAND_H
#define NORMALIZESHAREDCOMMAND_H

/*
 *  normalizesharedcommand.h
 *  Mothur
 *
 *  Created by westcott on 9/15/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */
 
#include "command.hpp"
#include "inputdata.h"
#include "sharedrabundvector.h"

class GlobalData;

class NormalizeSharedCommand : public Command {

public:
	NormalizeSharedCommand(string);
	NormalizeSharedCommand();
	~NormalizeSharedCommand();
	vector<string> getRequiredParameters();
	vector<string> getValidParameters();
	vector<string> getRequiredFiles();
	map<string, vector<string> > getOutputFiles() { return outputTypes; }
	int execute();
	void help();
	
private:
	GlobalData* globaldata;
	InputData* input;
	vector<SharedRAbundVector*> lookup;
	vector<SharedRAbundFloatVector*> lookupFloat;
	
	bool abort, allLines, pickedGroups, makeRelabund;
	set<string> labels; //holds labels to be used
	string groups, label, outputDir, method;
	int norm;
	vector<string> Groups, outputNames;
	map<string, vector<string> > outputTypes;
	
	int normalize(vector<SharedRAbundVector*>&, ofstream&);
	int normalize(vector<SharedRAbundFloatVector*>&, ofstream&);
	int eliminateZeroOTUS(vector<SharedRAbundVector*>&);
	int eliminateZeroOTUS(vector<SharedRAbundFloatVector*>&);

};

#endif

