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

class NormalizeSharedCommand : public Command {

public:
	NormalizeSharedCommand(string);
	NormalizeSharedCommand();
	~NormalizeSharedCommand() {}
	
	vector<string> setParameters();
	string getCommandName()			{ return "normalize.shared";		}
	string getCommandCategory()		{ return "OTU-Based Approaches";	}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "http://www.mothur.org/wiki/Normalize.shared"; }
	string getDescription()		{ return "normalize samples in a shared or relabund file"; }

	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
	
private:
	InputData* input;
	vector<SharedRAbundVector*> lookup;
	vector<SharedRAbundFloatVector*> lookupFloat;
	
	bool abort, allLines, pickedGroups, makeRelabund;
	set<string> labels; //holds labels to be used
	string groups, label, outputDir, method, sharedfile, relabundfile, format, inputfile;
	int norm;
	vector<string> Groups, outputNames;
		
	int normalize(vector<SharedRAbundVector*>&);
	int normalize(vector<SharedRAbundFloatVector*>&);
	int eliminateZeroOTUS(vector<SharedRAbundVector*>&);
	int eliminateZeroOTUS(vector<SharedRAbundFloatVector*>&);

};

#endif

