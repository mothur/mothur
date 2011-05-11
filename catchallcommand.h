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
	string getCitation() { return "Bunge, J. (2011). Estimating the number of species with CatchAll. Forthcoming in Proceedings of the Pacific Symposium on Biocomputing 2011.\nhttp://www.northeastern.edu/catchall/index.html http://www.mothur.org/wiki/Catchall"; }

	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	string outputDir, sharedfile, sabundfile, format, path, savedOutputDir;
	bool abort, allLines;
	set<string> labels;
	vector<string> outputNames;
	vector<string> groups;
	
	string process(SAbundVector*, string);
	int createSummaryFile(string, string, ofstream&); 
	vector<string> parseSharedFile(string);
	string combineSummmary(vector<string>&);
};

/****************************************************************************/

#endif


