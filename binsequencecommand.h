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

class BinSeqCommand : public Command {
	
public:
	BinSeqCommand(string);	
	BinSeqCommand();
	~BinSeqCommand();
	
	vector<string> setParameters();
	string getCommandName()			{ return "bin.seqs";			}
	string getCommandCategory()		{ return "Sequence Processing"; }
	string getOutputFileNameTag(string, string);
	string getHelpString();	
	string getCitation() { return "http://www.mothur.org/wiki/Bin.seqs"; }
	string getDescription()		{ return "maps sequences to otus"; }
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }		
	
private:
	
	ListVector* list;
	InputData* input;
	FastaMap* fasta;
	GroupMap* groupMap;
	bool abort, allLines;
	set<string> labels; //holds labels to be used
	string filename, fastafile, listfile, namesfile, groupfile, label, outputDir;
	ofstream out;
	ifstream in, inNames;
	vector<string> outputNames;
	
	void readNamesFile();
	int process(ListVector*);
};

#endif
