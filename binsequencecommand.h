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
#include "readotu.h"
#include "fastamap.h"
#include "groupmap.h"


class GlobalData;

class BinSeqCommand : public Command {
	
public:
	BinSeqCommand();	
	~BinSeqCommand();
	int execute();	
	
private:
	GlobalData* globaldata;
	ListVector* list;
	ReadOTUFile* read;
	InputData* input;
	FastaMap* fasta;
	GroupMap* groupMap;
	string filename, fastafile, namesfile, groupfile;
	ofstream out;
	ifstream in, inNames;
	
	void readNamesFile();
	int process(ListVector*, int);
};

#endif
