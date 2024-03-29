#ifndef CLUSTERFRAGMENTSCOMMAND_H
#define CLUSTERFRAGMENTSCOMMAND_H

/*
 *  clusterfragmentscommand.h
 *  Mothur
 *
 *  Created by westcott on 9/23/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */


#include "command.hpp"
#include "sequence.hpp"
#include "counttable.h"

/************************************************************/
struct seqRNode {
	int numIdentical;
	int length;
	Sequence seq;
	string names;
	bool active;
    seqRNode() = default;
	seqRNode(int n, Sequence s, string nm, int l) : numIdentical(n), seq(s), names(nm), active(1), length(l) {}
	~seqRNode() = default;
};
/************************************************************/

class ClusterFragmentsCommand : public Command {
	
public:
	ClusterFragmentsCommand(string);
	~ClusterFragmentsCommand() = default;
	
	vector<string> setParameters();
	string getCommandName()			{ return "cluster.fragments";		}
	string getCommandCategory()		{ return "Sequence Processing";		}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "http://www.mothur.org/wiki/Cluster.fragments"; }
	string getDescription()		{ return "creates a namesfile with sequences that are a fragment of a larger sequence"; }
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
    CountTable ct;
	bool abort;
	string fastafile, namefile, countfile;
	int diffs, percent;
	vector<seqRNode> alignSeqs; 
	map<string, string> names; //represents the names file first column maps to second column
	map<string, int> sizes;  //this map a seq name to the number of identical seqs in the names file
	map<string, int>::iterator itSize; 
	vector<string> outputNames;
	
	int readFASTA();
	void readNameFile();
	void printData(string, string); //fasta filename, names file name
	bool isFragment(string, string);
	
};

/************************************************************/

#endif

