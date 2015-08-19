#ifndef DECONVOLUTECOMMAND_H
#define DECONVOLUTECOMMAND_H
/*
 *  deconvolute.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 1/21/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "command.hpp"
#include "fastamap.h"
#include "counttable.h"

/* The unique.seqs command reads a fasta file, finds the duplicate sequences and outputs a names file
	containing 2 columns.  The first being the groupname and the second the list of identical sequence names. */ 


class DeconvoluteCommand : public Command {

public:
	DeconvoluteCommand(string);
	DeconvoluteCommand();
	~DeconvoluteCommand() {}
	
	vector<string> setParameters();
	string getCommandName()			{ return "unique.seqs";		}
	string getCommandCategory()		{ return "Sequence Processing";		}
    
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "http://www.mothur.org/wiki/Unique.seqs"; }
	string getDescription()		{ return "creates a fasta containing the unique sequences as well as a namesfile with the names each sequence represents"; }

	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
	
private:
	string fastafile, namefile, outputDir, countfile, format;
	vector<string> outputNames;

	bool abort;
};

#endif
