#ifndef PCACOMMAND_H
#define PCACOMMAND_H

/*
 *  pcacommand.h
 *  mothur
 *
 *  Created by westcott on 1/7/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "command.hpp"
#include "linearalgebra.h"
#include "globaldata.hpp"


/*****************************************************************/
class PCACommand : public Command {
	
public:
	PCACommand(string);	
	PCACommand();
	~PCACommand();
	vector<string> getRequiredParameters();
	vector<string> getValidParameters();
	vector<string> getRequiredFiles();
	map<string, vector<string> > getOutputFiles() { return outputTypes; }
	int execute();	
	void help();
	
private:
	GlobalData* globaldata;
	
	bool abort, metric;
	string outputDir, mode, inputFile, label, groups;
	vector<string> outputNames, Groups;
	set<string> labels;
	map<string, vector<string> > outputTypes;
	LinearAlgebra linearCalc;
	
	//vector< vector<double> > createMatrix(vector<SharedRAbundFloatVector*>);
	int process(vector<SharedRAbundFloatVector*>&);
	void output(string, vector<string>, vector<vector<double> >&, vector<double>);
	
};

/*****************************************************************/

#endif


