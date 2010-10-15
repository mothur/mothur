#ifndef NOCOMMAND_H
#define NOCOMMAND_H
/*
 *  nocommand.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

/* This command is run if the user enters an invalid command. */

#include "command.hpp"
#include "commandfactory.hpp"

class NoCommand : public Command {

public:
	NoCommand(string);
	NoCommand() {}
	~NoCommand();
	vector<string> getRequiredParameters();
	vector<string> getValidParameters();
	vector<string> getRequiredFiles();
	map<string, vector<string> > getOutputFiles() { return outputTypes; }
	int execute();
	void help() {}
	
private:
	vector<string> outputNames;
	map<string, vector<string> > outputTypes;
		
};

#endif
