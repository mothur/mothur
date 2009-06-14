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
	~LibShuffCommand(){};
	int execute();
	void help();	

private:
	vector<string> groupNames;
	
	void setGroups();
	void printCoverageFile();
	void printSummaryFile();

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
	string outputFile, groups, userform, savegroups;
	vector<string> Groups; //holds groups to be used
};

#endif
