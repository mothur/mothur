#ifndef REVERSECOMMAND_H
#define REVERSECOMMAND_H

/*
 *  reversecommand.h
 *  Mothur
 *
 *  Created by Pat Schloss on 6/6/09.
 *  Copyright 2009 Patrick D. Schloss. All rights reserved.
 *
 */

#include "command.hpp"

class ReverseSeqsCommand : public Command {
public:
	ReverseSeqsCommand(string);
	~ReverseSeqsCommand() {}
	
	vector<string> setParameters();
	string getCommandName()			{ return "reverse.seqs";			}
	string getCommandCategory()		{ return "Sequence Processing";		}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "http://www.mothur.org/wiki/Reverse.seqs"; }
	string getDescription()		{ return "outputs a fasta file containing the reverse-complements"; }

	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:

	bool abort;
	string fastaFileName, qualFileName;
	vector<string> outputNames;
	
};

#endif
