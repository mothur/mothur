#ifndef GETCURRENTCOMMAND_H
#define GETCURRENTCOMMAND_H

/*
 *  getcurrentcommand.h
 *  Mothur
 *
 *  Created by westcott on 3/16/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "command.hpp"

class GetCurrentCommand : public Command {

	public:
		GetCurrentCommand(string);
		GetCurrentCommand();
		~GetCurrentCommand();
		vector<string> getRequiredParameters();
		vector<string> getValidParameters();
		vector<string> getRequiredFiles();
		map<string, vector<string> > getOutputFiles() { return outputTypes; }
		int execute();
		void help();
		
	private:
		
		vector<string> outputNames;
		map<string, vector<string> > outputTypes;
		bool abort;
	
		string clearTypes;
		vector<string> types;
		
};

#endif

