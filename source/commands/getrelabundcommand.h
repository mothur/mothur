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
#include "sharedrabundvector.h"


class GetRelAbundCommand : public Command {

public:
	GetRelAbundCommand(string);
	GetRelAbundCommand();
	~GetRelAbundCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "get.relabund";			}
	string getCommandCategory()		{ return "OTU-Based Approaches";	}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "http://www.mothur.org/wiki/Get.relabund"; }
	string getDescription()		{ return "calculates the relative abundance of each OTU in a sample"; }

	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	InputData* input;
	vector<SharedRAbundVector*> lookup;
	
	bool abort, allLines, pickedGroups;
	set<string> labels; //holds labels to be used
	string groups, label, outputDir, scale, sharedfile;
	vector<string> Groups, outputNames;
	
	int getRelAbundance(vector<SharedRAbundVector*>&, ofstream&);

};

#endif

