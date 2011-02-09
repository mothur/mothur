#ifndef MANTELCOMMAND_H
#define MANTELCOMMAND_H

/*
 *  mantelcommand.h
 *  mothur
 *
 *  Created by westcott on 2/9/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "command.hpp"
#include "linearalgebra.h"

class MantelCommand : public Command {
public:
	MantelCommand(string);
	MantelCommand();
	~MantelCommand();
	vector<string> getRequiredParameters();
	vector<string> getValidParameters();
	vector<string> getRequiredFiles();
	map<string, vector<string> > getOutputFiles() { return outputTypes; }
	int execute();
	void help();
	
private:
	
	string phylipfile1, phylipfile2, outputDir, method;
	bool abort;
	int iters;
	
	vector<string> outputNames;
	map<string, vector<string> > outputTypes;
};


#endif



