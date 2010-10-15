#ifndef SYSTEMCOMMAND_H
#define SYSTEMCOMMAND_H

/*
 *  systemcommand.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 7/8/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "command.hpp"


class SystemCommand : public Command {
	
	public:
	
		SystemCommand(string);	
		SystemCommand() {}
		~SystemCommand(){}
		vector<string> getRequiredParameters();
		vector<string> getValidParameters();
		vector<string> getRequiredFiles();
		map<string, vector<string> > getOutputFiles() { return outputTypes; }
		int execute();
		void help();	
		
	private:
		string command;
		bool abort;
		vector<string> outputNames;
		map<string, vector<string> > outputTypes;
		
};

#endif

