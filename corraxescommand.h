#ifndef CORRAXESCOMMAND_H
#define CORRAXESCOMMAND_H

/*
 *  corraxescommand.h
 *  Mothur
 *
 *  Created by westcott on 12/22/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "command.hpp"
#include "sharedrabundvector.h"
#include "sharedrabundfloatvector.h"
#include "inputdata.h"

class CorrAxesCommand : public Command {
public:
	CorrAxesCommand(string);
	CorrAxesCommand();
	~CorrAxesCommand();
	vector<string> getRequiredParameters();
	vector<string> getValidParameters();
	vector<string> getRequiredFiles();
	map<string, vector<string> > getOutputFiles() { return outputTypes; }
	int execute();
	void help();
	
private:
	GlobalData* globaldata;
	string axesfile, sharedfile, relabundfile, metadatafile, groups, label, inputFileName, outputDir, method;
	bool abort, pickedGroups;
	int numaxes;
	set<string> names;
	
	vector<string> outputNames, Groups;
	map<string, vector<string> > outputTypes;
	vector<SharedRAbundVector*> lookup;
	vector<SharedRAbundFloatVector*> lookupFloat;
	
	int getShared();
	int getSharedFloat();
	int convertToRelabund();
	int eliminateZeroOTUS(vector<SharedRAbundFloatVector*>&);
	map<string, vector<float> > readAxes();
	int calcPearson(map<string, vector<float> >&, ofstream&);
	int calcSpearman(map<string, vector<float> >&, ofstream&);
	
};


#endif


