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
	bool active;
	seqPNode() {}
	seqPNode(int n, Sequence s, string nm) : numIdentical(n), seq(s), names(nm), active(1) {}
	~seqPNode() {}
};
/************************************************************/

class PreClusterCommand : public Command {
	
public:
	PreClusterCommand(string);
	PreClusterCommand();
	~PreClusterCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "pre.cluster";				}
	string getCommandCategory()		{ return "Sequence Processing";		}
	string getHelpString();	
	string getCitation() { return "http://www.mothur.org/wiki/Pre.cluster"; }
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	int diffs, length;
	bool abort;
	string fastafile, namefile, outputDir;
	vector<seqPNode> alignSeqs; //maps the number of identical seqs to a sequence
	map<string, string> names; //represents the names file first column maps to second column
	map<string, int> sizes;  //this map a seq name to the number of identical seqs in the names file
	map<string, int>::iterator itSize; 
//	map<string, bool> active; //maps sequence name to whether it has already been merged or not.
	vector<string> outputNames;
	map<string, vector<string> > outputTypes;
	
	int readFASTA();
	void readNameFile();
	//int readNamesFASTA();
	int calcMisMatches(string, string);
	void printData(string, string); //fasta filename, names file name
};

/************************************************************/





#endif


