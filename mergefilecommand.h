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
	~MergeFileCommand();
	int execute();
	void help();	

private:
	vector<string> fileNames;
	string outputFileName;
	int numInputFiles;
	bool abort;
};

#endif
