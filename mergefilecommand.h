#ifndef MERGEFILECOMMAND_H
#define MERGEFILECOMMAND_H

/*
 *  mergefilecommand.h
 *  Mothur
 *
 *  Created by Pat Schloss on 6/14/09.
 *  Copyright 2009 Patrick D. Schloss. All rights reserved.
 *
 */

#include "mothur.h"
#include "command.hpp"

class MergeFileCommand : public Command {
public:
	MergeFileCommand(string);
	MergeFileCommand();
	~MergeFileCommand();
	vector<string> getRequiredParameters();
	vector<string> getValidParameters();
	vector<string> getRequiredFiles();
	map<string, vector<string> > getOutputFiles() { return outputTypes; }
	int execute();
	void help();	

private:
	vector<string> fileNames, outputNames;
	map<string, vector<string> > outputTypes;
	string outputFileName;
	int numInputFiles;
	bool abort;
};

#endif
