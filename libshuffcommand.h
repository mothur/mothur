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
#include "groupmap.h"


class LibShuffCommand : public Command {
	
public:
	LibShuffCommand(string);
	LibShuffCommand();	
	~LibShuffCommand(){};
	
	vector<string> setParameters();
	string getCommandName()			{ return "libshuff";				}
	string getCommandCategory()		{ return "Hypothesis Testing";		}
	string getHelpString();	
	
	int execute();
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	vector<string> groupNames;
	
	void setGroups();
	int printCoverageFile();
	int printSummaryFile();
	
	GroupMap* groupMap;
	FullMatrix* matrix;
	Libshuff* form;
	float cutOff, step;
	int numGroups, numComp, iters;
	string coverageFile, summaryFile, phylipfile, groupfile;
	vector<vector<int> > pValueCounts;
	vector<vector<double> > savedDXYValues;
	vector<vector<vector<double> > > savedMinValues;

	bool abort, sim;
	string outputFile, groups, userform, savegroups, outputDir;
	vector<string> Groups, outputNames; //holds groups to be used
};

#endif
