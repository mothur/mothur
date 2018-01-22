#ifndef SUMMARYQUALCOMMAND_H
#define SUMMARYQUALCOMMAND_H

/*
 *  summaryqualcommand.h
 *  Mothur
 *
 *  Created by westcott on 11/28/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */


#include "command.hpp"
#include "qualityscores.h"

/**************************************************************************************************/

class SummaryQualCommand : public Command {
public:
	SummaryQualCommand(string);
	SummaryQualCommand();
	~SummaryQualCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "summary.qual";			}
	string getCommandCategory()		{ return "Sequence Processing";		}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "http://www.mothur.org/wiki/Summary.qual"; }
	string getDescription()		{ return "summarize the quality of a set of sequences"; }
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	bool abort, hasNameMap;
	string qualfile, outputDir, namefile, countfile;
	vector<string> outputNames;
	map<string, int> nameMap;
	int processors;

	long long createProcessesCreateSummary(vector<int>&, vector<int>&, vector< vector<int> >&, string);
	int printQual(string, vector<int>&, vector<int>&, vector< vector<int> >&);
};
/**************************************************************************************************/

#endif

