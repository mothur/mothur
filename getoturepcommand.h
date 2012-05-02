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
#include "listvector.hpp"
#include "inputdata.h"
#include "fastamap.h"
#include "groupmap.h"
#include "readmatrix.hpp"
#include "formatmatrix.h"

typedef list<PCell>::iterator MatData;
typedef map<int, float> SeqMap;

struct repStruct {
		string name;
		int	bin;
		int size;
		string group;
		
		repStruct(){}
		repStruct(string n, int b, int s, string g) : name(n), bin(b), size(s), group(g) {}
		~repStruct() {}
};

class GetOTURepCommand : public Command {

public:
	GetOTURepCommand(string);
	GetOTURepCommand();
	~GetOTURepCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "get.oturep";				}
	string getCommandCategory()		{ return "OTU-Based Approaches";	}
	string getHelpString();	
	string getCitation() { return "http://www.mothur.org/wiki/Get.oturep"; }
	string getDescription()		{ return "gets a representative sequence for each OTU"; }

	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
	
private:
	ListVector* list;
	InputData* input;
	FastaMap* fasta;
	GroupMap* groupMap;
	ReadMatrix* readMatrix;
	FormatMatrix* formatMatrix;
	NameAssignment* nameMap;
	string filename, fastafile, listfile, namefile, groupfile, label, sorted, phylipfile, columnfile, distFile, format, outputDir, groups;
	ofstream out;
	ifstream in, inNames, inRow;
	bool abort, allLines, groupError, large, weighted;
	set<string> labels; //holds labels to be used
	map<string, int> nameToIndex;  //maps sequence name to index in sparsematrix
	map<string, string> nameFileMap;
	vector<string> outputNames, Groups;
	map<string, string> outputNameFiles;
	float cutoff;
	int precision;
	vector<SeqMap> seqVec;			// contains maps with sequence index and distance
									// for all distances related to a certain sequence
	vector<int> rowPositions;

	void readNamesFile();
	void readNamesFile(bool);
	int process(ListVector*);
	SeqMap getMap(int);
	string findRep(vector<string>); 	// returns the name of the "representative" sequence of given bin or subset of a bin, for groups
	int processNames(string, string);
	int processFastaNames(string, string);
};

#endif

