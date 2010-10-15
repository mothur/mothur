#ifndef LIBSHUFFCOMMAND_H
#define LIBSHUFFCOMMAND_H

/*
 *  libshuffcommand.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 3/9/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "command.hpp"
#include "fullmatrix.h"
#include "libshuff.h"


class GlobalData;

class LibShuffCommand : public Command {
	
public:
	LibShuffCommand(string);
	LibShuffCommand();	
	~LibShuffCommand(){};
	vector<string> getRequiredParameters();
	vector<string> getValidParameters();
	vector<string> getRequiredFiles();
	map<string, vector<string> > getOutputFiles() { return outputTypes; }
	int execute();
	void help();	

private:
	vector<string> groupNames;
	
	void setGroups();
	int printCoverageFile();
	int printSummaryFile();

	GlobalData* globaldata;
	FullMatrix* matrix;
	Libshuff* form;
	float cutOff, step;
	int numGroups, numComp, iters;
	string coverageFile, summaryFile;
	vector<vector<int> > pValueCounts;
	vector<vector<double> > savedDXYValues;
	vector<vector<vector<double> > > savedMinValues;

	bool abort;
	string outputFile, groups, userform, savegroups, outputDir;
	vector<string> Groups, outputNames; //holds groups to be used
	map<string, vector<string> > outputTypes;
};

#endif
