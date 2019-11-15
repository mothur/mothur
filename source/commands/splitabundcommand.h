#ifndef SPLITABUNDCOMMAND_H
#define SPLITABUNDCOMMAND_H

/*
 *  splitabundcommand.h
 *  Mothur
 *
 *  Created by westcott on 5/17/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */


/* split.abund - given a list or name file and a number (cutoff), make two files - *rare* and *abund* 
- where rare has data for otus that have fewer sequences than the cutoff and abund has data for otus 
that have as many or more sequences as the cutoff. 
also allow an option where a user can give a group file with the list or names file and split the group file into rare and abund. */


#include "command.hpp"
#include "groupmap.h"
#include "inputdata.h"
#include "listvector.hpp"
#include "sequence.hpp"
#include "counttable.h"

/***************************************************************************************/

class SplitAbundCommand : public Command {
	
public:
	SplitAbundCommand(string);	
	SplitAbundCommand();
	~SplitAbundCommand();
	
	vector<string> setParameters();
	string getCommandName()			{ return "split.abund";				}
	string getCommandCategory()		{ return "OTU-Based Approaches";	}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "http://www.mothur.org/wiki/Split.abund"; }
	string getDescription()		{ return "split a list, name, group or fasta file based on abundance"; }
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	int splitList();
    int splitCount();
    int splitNames();
    int process(ListVector*);
	int writeList(ListVector*, string, int); 
	vector<string> writeAccnos(string, set<string>, set<string>);
    
	vector<string> outputNames;
    CountTable ct;
	
	string outputDir, listfile, namefile, groupfile, countfile, label, groups, fastafile, inputFile;
	set<string> labels;
	bool abort, allLines, accnos;
	float cutoff;
};

/***************************************************************************************/

#endif


