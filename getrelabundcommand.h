#ifndef GETRELABUNDCOMMAND_H
#define GETRELABUNDCOMMAND_H

/*
 *  getrelabundcommand.h
 *  Mothur
 *
 *  Created by westcott on 6/21/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */
 
#include "command.hpp"
#include "inputdata.h"
#include "readotu.h"
#include "sharedrabundvector.h"

class GlobalData;

class GetRelAbundCommand : public Command {

public:
	GetRelAbundCommand(string);
	GetRelAbundCommand();
	~GetRelAbundCommand();
	vector<string> getRequiredParameters();
	vector<string> getValidParameters();
	vector<string> getRequiredFiles();
	map<string, vector<string> > getOutputFiles() { return outputTypes; }
	int execute();
	void help();
	
private:
	GlobalData* globaldata;
	ReadOTUFile* read;
	InputData* input;
	vector<SharedRAbundVector*> lookup;
	
	bool abort, allLines, pickedGroups;
	set<string> labels; //holds labels to be used
	string groups, label, outputDir, scale;
	vector<string> Groups, outputNames;
	map<string, vector<string> > outputTypes;
	
	int getRelAbundance(vector<SharedRAbundVector*>&, ofstream&);
	int eliminateZeroOTUS(vector<SharedRAbundVector*>& thislookup);

};

#endif

