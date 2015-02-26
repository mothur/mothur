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
#include "sharedrabundfloatvector.h"

/*****************************************************************/
class PCACommand : public Command {
	
public:
	PCACommand(string);	
	PCACommand();
	~PCACommand() {}
	
	vector<string> setParameters();
	string getCommandName()			{ return "pca";					}
	string getCommandCategory()		{ return "Hypothesis Testing";	}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "McCune B, Grace JB, Urban DL (2002). Analysis of ecological communities. MjM Software Design: Gleneden Beach, OR. \nLegendre P, Legendre L (1998). Numerical Ecology. Elsevier: New York. \nhttp://www.mothur.org/wiki/Pca"; }
	string getDescription()		{ return "pca"; }

	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
	
private:

	bool abort, metric;
	string outputDir, mode, inputFile, label, groups, sharedfile, relabundfile;
	vector<string> outputNames, Groups;
	set<string> labels;
	LinearAlgebra linearCalc;
	
	//vector< vector<double> > createMatrix(vector<SharedRAbundFloatVector*>);
	int process(vector<SharedRAbundFloatVector*>&);
	void output(string, string, vector<string>, vector<vector<double> >&, vector<double>);
	
};

/*****************************************************************/

#endif


