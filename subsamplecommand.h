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
#include "globaldata.hpp"
#include "sharedrabundvector.h"
#include "listvector.hpp"
#include "rabundvector.hpp"
#include "inputdata.h"


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
	
	bool abort, pickedGroups, allLines;
	string listfile, groupfile, sharedfile, rabundfile, sabundfile, fastafile, namefile;
	set<string> labels; //holds labels to be used
	string groups, label, outputDir;
	vector<string> Groups, outputNames;
	map<string, vector<string> > outputTypes;
	int size;
	vector<string> names;
	map<string, vector<string> > nameMap;
	
	int eliminateZeroOTUS(vector<SharedRAbundVector*>&);
	int getSubSampleShared();
	int getSubSampleList();
	int getSubSampleRabund();
	int getSubSampleSabund();
	int getSubSampleFasta();
	int processShared(vector<SharedRAbundVector*>&, ofstream&);
	int processRabund(RAbundVector*&, ofstream&);
	int processSabund(SAbundVector*&, ofstream&);
	int processList(ListVector*&, ofstream&, set<string>&);
	int getNames();
	int readNames();
	
};

#endif

