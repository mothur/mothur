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

/***************************************************************************************/

class SplitGroupCommand : public Command {
	
public:
	SplitGroupCommand(string);	
	~SplitGroupCommand();
	int execute();	
	void help();

	
private:
	int readNames(); 
	int splitFasta(); 
	
	vector<string> outputNames;
	map<string, vector<string> > nameMap;
	map<string, vector<string> >::iterator itNames;
	GroupMap* groupMap;
	
	string outputDir, namefile, groupfile, groups, fastafile;
	vector<string> Groups;
	bool abort;
};

/***************************************************************************************/

#endif



