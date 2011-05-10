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
#include "sharedrabundfloatvector.h"
#include "inputdata.h"


class CorrAxesCommand : public Command {
public:
	CorrAxesCommand(string);
	CorrAxesCommand();
	~CorrAxesCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "corr.axes";				}
	string getCommandCategory()		{ return "Hypothesis Testing";		}
	string getHelpString();	
	string getCitation() { return "http://www.mothur.org/wiki/Corr.axes"; }
	
	int execute();
	void help() { m->mothurOut(getHelpString()); }	
private:

	string axesfile, sharedfile, relabundfile, metadatafile, groups, label, inputFileName, outputDir, method;
	bool abort, pickedGroups;
	int numaxes;
	set<string> names;
	
	vector<string> outputNames, Groups;
	vector<SharedRAbundFloatVector*> lookupFloat;
	vector<string> metadataLabels;
	
	int getSharedFloat(InputData*);
	int getMetadata();
	int eliminateZeroOTUS(vector<SharedRAbundFloatVector*>&);
	map<string, vector<float> > readAxes();
	int calcPearson(map<string, vector<float> >&, ofstream&);
	int calcSpearman(map<string, vector<float> >&, ofstream&);
	int calcKendall(map<string, vector<float> >&, ofstream&);
	
};


#endif


