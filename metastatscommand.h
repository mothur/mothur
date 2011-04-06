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
#include "sharedrabundvector.h"

class MetaStatsCommand : public Command {

public:
	MetaStatsCommand(string);
	MetaStatsCommand();
	~MetaStatsCommand() {}
	
	vector<string> setParameters();
	string getCommandName()			{ return "metastats";				}
	string getCommandCategory()		{ return "OTU-Based Approaches";	}
	string getHelpString();	
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	struct linePair {
		int start;
		int num;
		linePair(int i, int j) : start(i), num(j) {}
	};
	vector<linePair> lines;
	
	GroupMap* designMap;
	InputData* input;
	vector<SharedRAbundVector*> lookup;
		
	bool abort, allLines, pickedGroups;
	set<string> labels; //holds labels to be used
	string groups, label, outputDir, inputDir, designfile, sets, sharedfile;
	vector<string> Groups, outputNames, Sets;
	vector< vector<string> > namesOfGroupCombos;
	int iters, processors;
	float threshold;
	
	int process(vector<SharedRAbundVector*>&);
	int driver(int, int, vector<SharedRAbundVector*>&);
};

#endif

