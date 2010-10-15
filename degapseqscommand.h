#ifndef DEGAPSEQSCOMMAND_H
#define DEGAPSEQSCOMMAND_H

/*
 *  degapseqscommand.h
 *  Mothur
 *
 *  Created by westcott on 6/21/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */


#include "command.hpp"

class DegapSeqsCommand : public Command {
public:
	DegapSeqsCommand(string);
	DegapSeqsCommand();
	~DegapSeqsCommand();
	vector<string> getRequiredParameters();
	vector<string> getValidParameters();
	vector<string> getRequiredFiles();
	map<string, vector<string> > getOutputFiles() { return outputTypes; }
	int execute();
	void help();
	
private:

	bool abort;
	string fastafile, outputDir;
	vector<string> outputNames;
	vector<string> fastaFileNames;
	map<string, vector<string> > outputTypes;
};

#endif

