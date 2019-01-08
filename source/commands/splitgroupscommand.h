#ifndef SPLITGROUPSCOMMAND_H
#define SPLITGROUPSCOMMAND_H

/*
 *  splitgroupscommand.h
 *  Mothur
 *
 *  Created by westcott on 9/20/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */


/* split.groups - given a group file, split sequences and names files in to separate files *.group1.fasta and .group1.names. */


#include "command.hpp"
#include "groupmap.h"
#include "sequence.hpp"
#include "fastqread.h"

/***************************************************************************************/

class SplitGroupCommand : public Command {
	
public:
	SplitGroupCommand(string);	
	SplitGroupCommand();
	~SplitGroupCommand() {}
	
	vector<string> setParameters();
	string getCommandName()			{ return "split.groups";				}
	string getCommandCategory()		{ return "Sequence Processing";		}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "http://www.mothur.org/wiki/Split.group"; }
	string getDescription()		{ return "split a name or fasta file by group"; }

	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	vector<string> outputNames;
		
	string outputDir, namefile, groupfile, countfile, groups, fastafile, flowfile, fastqfile, format;
	vector<string> Groups;
	bool abort;
    
    int runNameGroup();
    int runCount();
    int splitFlow();
    int splitFastq();
};

/***************************************************************************************/

#endif



