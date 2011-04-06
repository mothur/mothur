#ifndef MATRIXOUTPUTCOMMAND_H
#define MATRIXOUTPUTCOMMAND_H

/*
 *  matrixoutputcommand.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 5/20/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */ 
#include "command.hpp"
#include "inputdata.h"
#include "groupmap.h"
#include "validcalculator.h"

// aka. dist.shared()

/* This command create a tree file for each similarity calculator at distance level, using various calculators to find the similiarity between groups. 
	The user can select the labels they wish to use as well as the groups they would like included.
	They can also use as many or as few calculators as they wish. */
	

class MatrixOutputCommand : public Command {
	
public:
	MatrixOutputCommand(string);
	MatrixOutputCommand();	
	~MatrixOutputCommand();
	
	vector<string> setParameters();
	string getCommandName()			{ return "dist.shared";				}
	string getCommandCategory()		{ return "OTU-Based Approaches";	}
	string getHelpString();	
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	void printSims(ostream&);
	int process(vector<SharedRAbundVector*>);
	
	vector<Calculator*> matrixCalculators;
	vector< vector<float> > simMatrix;
	InputData* input;
	vector<SharedRAbundVector*> lookup;
	string exportFileName, output, sharedfile;
	int numGroups;
	ofstream out;

	bool abort, allLines;
	set<string> labels; //holds labels to be used
	string outputFile, calc, groups, label, outputDir;
	vector<string>  Estimators, Groups, outputNames; //holds estimators to be used
};
	
	
#endif

