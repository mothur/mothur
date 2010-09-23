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

/************************************************************/
struct seqRNode {
	int numIdentical;
	int length;
	Sequence seq;
	string names;
	bool active;
	seqRNode() {}
	seqRNode(int n, Sequence s, string nm, int l) : numIdentical(n), seq(s), names(nm), active(1), length(l) {}
	~seqRNode() {}
};
/************************************************************/

class ClusterFragmentsCommand : public Command {
	
public:
	ClusterFragmentsCommand(string);
	~ClusterFragmentsCommand();
	int execute();	
	void help();
	
private:
	bool abort;
	string fastafile, namefile, outputDir;
	vector<seqRNode> alignSeqs; 
	map<string, string> names; //represents the names file first column maps to second column
	map<string, int> sizes;  //this map a seq name to the number of identical seqs in the names file
	map<string, int>::iterator itSize; 
	
	int readFASTA();
	void readNameFile();
	void printData(string, string); //fasta filename, names file name
	
};

/************************************************************/

#endif

