#ifndef CONSENSUSSEQSCOMMAND_H
#define CONSENSUSSEQSCOMMAND_H

/*
 *  consensusseqscommand.h
 *  Mothur
 *
 *  Created by westcott on 11/23/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */


#include "command.hpp"
#include "listvector.hpp"

class ConsensusSeqsCommand : public Command {
public:
	ConsensusSeqsCommand(string);
	ConsensusSeqsCommand();
	~ConsensusSeqsCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "consensus.seqs";		}
	string getCommandCategory()		{ return "Sequence Processing"; }
	string getHelpString();	
	string getCitation() { return "http://www.mothur.org/wiki/Consensus.seqs"; }
	string getDescription()		{ return "create a consensus sequence for each OTU or for a fasta file"; }

	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	
	bool abort, allLines;
	string fastafile, listfile, namefile, label, outputDir;
	set<string> labels;
	vector<string> outputNames;
	map<string, string> fastaMap;
	map<string, string> nameMap;
	map<string, string> nameFileMap;
	int cutoff;
	
	int readFasta();
	int readNames();
	int processList(ListVector*&);
	string getConsSeq(string, ofstream&, string&, int);
	char getBase(vector<int>, int);
};

#endif




