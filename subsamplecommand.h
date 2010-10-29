#ifndef SUBSAMPLECOMMAND_H
#define SUBSAMPLECOMMAND_H

/*
 *  subsamplecommand.h
 *  Mothur
 *
 *  Created by westcott on 10/27/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */
 
#include "command.hpp"
#include "inputdata.h"
#include "readotu.h"
#include "sharedrabundvector.h"

class GlobalData;

class SubSampleCommand : public Command {

public:
	SubSampleCommand(string);
	SubSampleCommand();
	~SubSampleCommand();
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
	ListVector* list;
	RAbundVector* rabund;
	SAbundVector* sabund;
	
	bool abort, allLines, pickedGroups;
	set<string> labels; //holds labels to be used
	string groups, label, outputDir;
	vector<string> Groups, outputNames;
	map<string, vector<string> > outputTypes;
	
	int eliminateZeroOTUS(vector<SharedRAbundVector*>&);
	int getSubSampleShared(vector<SharedRAbundVector*>&, ofstream&);
};

#endif

