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
#include "groupmap.h"

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
	GroupMap* groupMap;
	string filename, fastafile, namesfile, groupfile;
	ofstream out;
	ifstream in, inNames;
	bool groupError;
	
	 
	map<string, int> nameToIndex;  //maps sequence name to index in sparsematrix
	map<int, string>::iterator it;
	map<int, string>::iterator it2;
	map<string, int>::iterator it3;
	
	void readNamesFile();
	int process(ListVector*);
	string FindRep(int, string&, ListVector*); // returns name of "representative" sequence of given bin. //and fill a string containing the groups in that bin if a groupfile is given

};

#endif

