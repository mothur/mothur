#ifndef BINSEQCOMMAND_H
#define BINSEQCOMMAND_H
/*
 *  binsequencecommand.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 4/3/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

/* The bin.seqs command outputs a .fasta file for each distance you specify appending the OTU number to each name. */

#include "command.hpp"
#include "inputdata.h"
#include "listvector.hpp"
#include "fastamap.h"
#include "groupmap.h"
#include "counttable.h"
#include "getseqscommand.h"

class BinSeqCommand : public Command {
	
public:
	BinSeqCommand(string);	
	~BinSeqCommand();
	
	vector<string> setParameters();
	string getCommandName()			{ return "bin.seqs";			}
	string getCommandCategory()		{ return "Sequence Processing"; }
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "http://www.mothur.org/wiki/Bin.seqs"; }
	string getDescription()		{ return "maps sequences to otus"; }
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }		
	
private:
	CountTable ct;
	bool abort, allLines;
	set<string> labels; //holds labels to be used
	string filename, fastafile, listfile, countfile, label;
	vector<string> outputNames;
	
	int process(ListVector*, FastaMap&);
};

#endif
