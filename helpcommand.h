#ifndef HELPCOMMAND_H
#define HELPCOMMAND_H
/*
 *  helpcommand.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

/* This class is designed to aid the user in running mothur. */

#include "command.hpp"
#include "commandfactory.hpp"

class HelpCommand : public Command {
	
public:
	HelpCommand(string);
	HelpCommand() {}
	~HelpCommand();
	vector<string> getRequiredParameters();
	vector<string> getValidParameters();
	vector<string> getRequiredFiles();
	map<string, vector<string> > getOutputFiles() { return outputTypes; }
	int execute();
	void help() {};
	
private:
	CommandFactory* validCommands;
	vector<string> outputNames;
	map<string, vector<string> > outputTypes;
	
	
private:
		
};
 
#endif
