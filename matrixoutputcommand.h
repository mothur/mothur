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
#include "readotu.h"
#include "validcalculator.h"

/* This command create a tree file for each similarity calculator at distance level, using various calculators to find the similiarity between groups. 
	The user can select the lines or labels they wish to use as well as the groups they would like included.
	They can also use as many or as few calculators as they wish. */
	
class GlobalData;

class MatrixOutputCommand : public Command {
	
public:
	MatrixOutputCommand();	
	~MatrixOutputCommand();
	int execute();	
	
private:
	void printSims(ostream&);
	
	GlobalData* globaldata;
	ReadOTUFile* read;
	vector<Calculator*> matrixCalculators;
	vector< vector<float> > simMatrix;
	InputData* input;
	ValidCalculators* validCalculator;
	vector<SharedRAbundVector*> lookup;
	string exportFileName;
	int numGroups;
	ofstream out;

};
	
	
#endif

