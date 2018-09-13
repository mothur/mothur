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
#include "counttable.h"
#include "optimatrix.h"
#include "nameassignment.hpp"
#include "countseqscommand.h"

typedef map<int, float> SeqMap;

struct repStruct {
		string name;
		string bin;
        int simpleBin;
		int size;
		string group;
		
		repStruct(){}
		repStruct(string n, string b, int sb, int s, string g) : name(n), bin(b), size(s), group(g), simpleBin(sb) { }
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
    map<string, int> nameMap;
    OptiData* matrix;
    CountTable ct;
	string filename, fastafile, listfile, namefile, groupfile, label, sorted, phylipfile, countfile, columnfile, distFile, format, outputDir, groups, method;
	ofstream out;
	ifstream in, inNames, inRow;
	bool abort, allLines, groupError, weighted, hasGroups;
	set<string> labels; //holds labels to be used
	//map<string, string> nameFileMap;
	vector<string> outputNames, Groups;
	map<string, string> outputNameFiles;
	float cutoff;
	int precision;
    
	

	void readNamesFile(FastaMap*&);
	//void readNamesFile(bool);
	int process(ListVector*);
	string findRep(vector<string>, string); 	// returns the name of the "representative" sequence of given bin or subset of a bin, for groups
    string findRepAbund(vector<string>, string);
	int processNames(string, string);
	int processFastaNames(string, string, FastaMap*&);
    int readDist();
};

#endif

