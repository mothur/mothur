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
#include "readmatrix.hpp"
#include "fastamap.h"


class GlobalData;

class BinSeqCommand : public Command {
	
public:
	BinSeqCommand();	
	~BinSeqCommand();
	int execute();	
	
private:
	GlobalData* globaldata;
	ListVector* list;
	ReadMatrix* read;
	InputData* input;
	FastaMap* fasta;
	string filename, fastafile, namesfile;
	ofstream out;
	ifstream in, inNames;
	
	void readNamesFile();
};

#endif
