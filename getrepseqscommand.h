#ifndef GETREPSEQSCOMMAND_H
#define GETREPSEQSCOMMAND_H
/*
 *  getrepseqscommand.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 5/19/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */


#include "command.hpp"
#include "inputdata.h"
#include "listvector.hpp"
#include "readotu.h"
#include "fastamap.h"
#include "groupmap.h"


class GlobalData;

class GetRepSeqsCommand : public Command {
	
public:
	GetRepSeqsCommand();	
	~GetRepSeqsCommand();
	int execute();	
	
private:
	GlobalData* globaldata;
	ListVector* list;
	ReadOTUFile* read;
	GroupMap* groupMap;
	InputData* input;
	FastaMap* fasta;
	string filename, fastafile, namesfile;
	map<string, ofstream*> filehandles;
	map<string, ofstream*>::iterator it;
	map<string, bool> used;  //group, if it had any unique otus
	map<string, bool>::iterator it2;
	map<string, string> seq;
	map<string, string>::iterator it3;
	ifstream in, inNames;
	
	void readNamesFile();
	void removeFiles(string);
};

#endif
