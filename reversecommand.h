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
	ReverseSeqsCommand();
	~ReverseSeqsCommand() {}
	
	vector<string> setParameters();
	string getCommandName()			{ return "reverse.seqs";			}
	string getCommandCategory()		{ return "Sequence Processing";		}
	string getHelpString();	
	string getCitation() { return "http://www.mothur.org/wiki/Reverse.seqs"; }
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:

	bool abort;
	string fastaFileName, qualFileName, outputDir;
	vector<string> outputNames;
	
};

#endif
