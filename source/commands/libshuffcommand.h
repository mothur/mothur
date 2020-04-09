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
	~LibShuffCommand(){};
	
	vector<string> setParameters();
	string getCommandName()			{ return "libshuff";				}
	string getCommandCategory()		{ return "Hypothesis Testing";		}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "Singleton DR, Furlong MA, Rathbun SL, Whitman WB (2001). Quantitative comparisons of 16S rRNA gene sequence libraries from environmental samples. Appl Environ Microbiol 67: 4374-6. \nSchloss PD, Larget BR, Handelsman J (2004). Integration of microbial ecology and statistics: a test to compare gene libraries. Appl Environ Microbiol 70: 5485-92. \nhttp://www.mothur.org/wiki/Libshuff"; }
	string getDescription()		{ return "a generic test that describes whether two or more communities have the same structure using the Cramer-von Mises test statistic"; }

	int execute();
	void help() { m->mothurOut(getHelpString()); }	
	
private:
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
