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
#include "counttable.h"

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
    string getOutputPattern(string);	
	string getCitation() { return "http://www.mothur.org/wiki/Get.oturep"; }
	string getDescription()		{ return "gets a representative sequence for each OTU"; }

	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
	
private:
	ListVector* list;
	GroupMap* groupMap;
	ReadMatrix* readMatrix;
	FormatMatrix* formatMatrix;
	NameAssignment* nameMap;
    CountTable ct;
	string filename, fastafile, listfile, namefile, groupfile, label, sorted, phylipfile, countfile, columnfile, distFile, format, outputDir, groups, method;
	ofstream out;
	ifstream in, inNames, inRow;
	bool abort, allLines, groupError, large, weighted, hasGroups;
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

	void readNamesFile(FastaMap*&);
	void readNamesFile(bool);
	int process(ListVector*);
	SeqMap getMap(int);
	string findRep(vector<string>, string); 	// returns the name of the "representative" sequence of given bin or subset of a bin, for groups
    string findRepAbund(vector<string>, string);
	int processNames(string, string);
	int processFastaNames(string, string, FastaMap*&);
    int readDist();
};

#endif

