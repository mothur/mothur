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
#include "readotu.h"
#include "sharedrabundvector.h"

class GlobalData;

class NormalizeSharedCommand : public Command {

public:
	NormalizeSharedCommand(string);
	~NormalizeSharedCommand();
	int execute();
	void help();
	
private:
	GlobalData* globaldata;
	ReadOTUFile* read;
	InputData* input;
	vector<SharedRAbundVector*> lookup;
	
	bool abort, allLines, pickedGroups;
	set<string> labels; //holds labels to be used
	string groups, label, outputDir, method;
	int norm;
	vector<string> Groups;
	
	int normalize(vector<SharedRAbundVector*>&, ofstream&);
	int eliminateZeroOTUS(vector<SharedRAbundVector*>& thislookup);

};

#endif

