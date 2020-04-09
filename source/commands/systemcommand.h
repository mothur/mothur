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
		~SystemCommand(){}
	
		vector<string> setParameters();
		string getCommandName()			{ return "system";		}
		string getCommandCategory()		{ return "General";		}
        string getHelpString();	
        string getOutputPattern(string){ return ""; }	
		string getCitation() { return "http://www.mothur.org/wiki/System"; }
		string getDescription()		{ return "execute system commands from within mothur"; }

		int execute(); 
		void help() { m->mothurOut(getHelpString()); }	
	
	private:
		string command;
		bool abort;
		vector<string> outputNames;
		
		
};

#endif

