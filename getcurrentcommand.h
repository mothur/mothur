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
		~GetCurrentCommand() {}
	
		vector<string> setParameters();
		string getCommandName()			{ return "get.current";	}
		string getCommandCategory()		{ return "General";		}
		string getHelpString();	
	
		int execute(); 
		void help() { m->mothurOut(getHelpString()); }	
	
	
	private:
		
		vector<string> outputNames;
		bool abort;
	
		string clearTypes;
		vector<string> types;
		
};

#endif

