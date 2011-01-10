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
	MetaStatsCommand();
	~MetaStatsCommand();
	vector<string> getRequiredParameters();
	vector<string> getValidParameters();
	vector<string> getRequiredFiles();
	map<string, vector<string> > getOutputFiles() { return outputTypes; }
	int execute();
	void help();
	
private:
	struct linePair {
		int start;
		int num;
		linePair(int i, int j) : start(i), num(j) {}
	};
	vector<linePair> lines;
	
	GlobalData* globaldata;
	GroupMap* designMap;
	ReadOTUFile* read;
	InputData* input;
	vector<SharedRAbundVector*> lookup;
	map<string, vector<string> > outputTypes;
	
	bool abort, allLines, pickedGroups;
	set<string> labels; //holds labels to be used
	string groups, label, outputDir, inputDir, designfile, sets;
	vector<string> Groups, outputNames, Sets;
	vector< vector<string> > namesOfGroupCombos;
	int iters, processors;
	float threshold;
	
	int process(vector<SharedRAbundVector*>&);
	int driver(int, int, vector<SharedRAbundVector*>&);
};

#endif

