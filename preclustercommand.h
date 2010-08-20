#ifndef PRECLUSTERCOMMAND_H
#define PRECLUSTERCOMMAND_H


/*
 *  preclustercommand.h
 *  Mothur
 *
 *  Created by westcott on 12/21/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */


#include "command.hpp"
#include "sequence.hpp"

/************************************************************/
struct seqPNode {
	int numIdentical;
	Sequence seq;
	string names;
	seqPNode() {}
	seqPNode(int n, Sequence s, string nm) : numIdentical(n), seq(s), names(nm) {}
	~seqPNode() {}
};
/************************************************************/

class PreClusterCommand : public Command {
	
public:
	PreClusterCommand(string);
	~PreClusterCommand();
	int execute();	
	void help();
	
private:
	int diffs, length;
	bool abort;
	string fastafile, namefile, outputDir;
	list<seqPNode> alignSeqs; //maps the number of identical seqs to a sequence
	map<string, string> names; //represents the names file first column maps to second column
	map<string, int> sizes;  //this map a seq name to the number of identical seqs in the names file
	map<string, int>::iterator itSize; 
//	map<string, bool> active; //maps sequence name to whether it has already been merged or not.
	
	int readFASTA();
	void readNameFile();
	//int readNamesFASTA();
	int calcMisMatches(string, string);
	void printData(ofstream&, ofstream&, seqPNode); //fasta filename, names file name
};

/************************************************************/





#endif


