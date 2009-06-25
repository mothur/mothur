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

typedef list<PCell>::iterator MatData;
typedef map<int, float> SeqMap;

class GetOTURepCommand : public Command {

public:
	GetOTURepCommand(string);
	~GetOTURepCommand();
	int execute();
	void help();

private:
	GlobalData* globaldata;
	SparseMatrix* matrix;
	ListVector* list;
	ReadOTUFile* read;
	InputData* input;
	FastaMap* fasta;
	GroupMap* groupMap;
	string filename, fastafile, listfile, namesfile, groupfile, line, label;
	ofstream out;
	ifstream in, inNames;
	bool groupError;

	bool abort, allLines;
	set<int> lines; //hold lines to be used
	set<string> labels; //holds labels to be used
	map<string, int> nameToIndex;  //maps sequence name to index in sparsematrix

	vector<SeqMap> seqVec;			// contains maps with sequence index and distance
									// for all distances related to a certain sequence


	void readNamesFile();
	int process(ListVector*);
	string findRep(int, string&, ListVector*, int&); 	// returns the name of the "representative" sequence of given bin, 
									// fills a string containing the groups in that bin if a groupfile is given,
									// and returns the number of sequences in the given bin

};

#endif

