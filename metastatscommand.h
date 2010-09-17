#ifndef METASTATSCOMMAND_H
#define METASTATSCOMMAND_H

/*
 *  metastatscommand.h
 *  Mothur
 *
 *  Created by westcott on 9/16/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */
 
#include "command.hpp"
#include "inputdata.h"
#include "readotu.h"
#include "sharedrabundvector.h"

class GlobalData;

class MetaStatsCommand : public Command {

public:
	MetaStatsCommand(string);
	~MetaStatsCommand();
	int execute();
	void help();
	
private:
	GlobalData* globaldata;
	GroupMap* designMap;
	ReadOTUFile* read;
	InputData* input;
	vector<SharedRAbundVector*> lookup;
	
	bool abort, allLines, pickedGroups;
	set<string> labels; //holds labels to be used
	string groups, label, outputDir, inputDir, designfile;
	vector<string> Groups, outputNames;
	int iters, g;
	float threshold;
	
	int process(vector<SharedRAbundVector*>&);
	int eliminateZeroOTUS(vector<SharedRAbundVector*>&);
};

#endif

