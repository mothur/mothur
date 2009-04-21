#ifndef GETOTUREPCOMMAND_H
#define GETOTUREPCOMMAND_H
/*
 *  getoturepcommand.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 4/6/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */
 
 /* The get.oturep command outputs a .fastarep file for each distance you specify, selecting one OTU representative for each bin. */

#include "command.hpp"
#include "globaldata.hpp"
#include "sparsematrix.hpp"
#include "listvector.hpp"
#include "inputdata.h"
#include "readotu.h"
#include "fastamap.h"


class GlobalData;

typedef list<PCell>::iterator MatData;

class GetOTURepCommand : public Command {
	
public:
	GetOTURepCommand();	
	~GetOTURepCommand();
	int execute();	
	
private:
	GlobalData* globaldata;
	SparseMatrix* matrix;
	ListVector* list;
	ListVector* listOfNames;
	ReadOTUFile* read;
	InputData* input;
	FastaMap* fasta;
	string filename, fastafile, namesfile;
	ofstream out;
	ifstream in, inNames;
	
	 
	map<string, int> nameToIndex;  //maps sequence name to index in sparsematrix
	map<int, string>::iterator it;
	map<int, string>::iterator it2;
	map<string, int>::iterator it3;
	
	void readNamesFile();
	string FindRep(int); // returns name of "representative" sequence of given bin.

};

#endif

