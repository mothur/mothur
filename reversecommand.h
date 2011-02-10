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
	~ReverseSeqsCommand();
	vector<string> getRequiredParameters();
	vector<string> getValidParameters();
	vector<string> getRequiredFiles();
	map<string, vector<string> > getOutputFiles() { return outputTypes; }
	int execute();
	void help();
	
private:

	bool abort;
	string fastaFileName, qualFileName, outputDir;
	vector<string> outputNames;
	map<string, vector<string> > outputTypes;
	
};

#endif
