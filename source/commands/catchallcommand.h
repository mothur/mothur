#ifndef CATCHALLCOMMAND_H
#define CATCHALLCOMMAND_H

/*
 *  catchallcommand.h
 *  Mothur
 *
 *  Created by westcott on 5/11/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "command.hpp"
#include "inputdata.h"
#include "sabundvector.hpp"

/* 
 citation goes here
 */ 

/****************************************************************************/

class CatchAllCommand : public Command {

public:

	CatchAllCommand(string);
	CatchAllCommand();
	~CatchAllCommand() {}
	
	vector<string> setParameters();
	string getCommandName()			{ return "catchall";			}
	string getCommandCategory()		{ return "OTU-Based Approaches";	}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "Bunge J, Woodard L, Bohning D, Foster JA, Connolly S, Allen HK (2012). Estimating population diversity with CatchAll. Bioinformatics  28:1045.\nhttp://www.northeastern.edu/catchall/index.html\nhttp://www.mothur.org/wiki/Catchall"; }
	string getDescription()		{ return "estimate number of species"; }
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	string outputDir, sharedfile, sabundfile, format, path, savedOutputDir;
	bool abort, allLines;
	set<string> labels;
	vector<string> outputNames;
	vector<string> Groups;
	
	string process(SAbundVector*, string);
	int createSummaryFile(string, string, ofstream&); 
	vector<string> parseSharedFile(string);
	string combineSummmary(vector<string>&);
};

/****************************************************************************/

#endif


