#ifndef COUNTGROUPSCOMMAND_H
#define COUNTGROUPSCOMMAND_H

/*
 *  countgroupscommand.h
 *  Mothur
 *
 *  Created by westcott on 8/9/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "command.hpp"

class CountGroupsCommand : public Command {
	
public:
	
	CountGroupsCommand(string);
	~CountGroupsCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "count.groups";			}
	string getCommandCategory()		{ return "Sequence Processing";		}
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "http://www.mothur.org/wiki/Count.groups"; }
	string getDescription()		{ return "counts the number of sequences in each group"; }
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
	
private:
	string sharedfile, groupfile, countfile, outputDir, groups, accnosfile;
	bool abort;
	vector<string> Groups;
    vector<string> outputNames;
};

#endif
