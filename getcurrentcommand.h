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
#include "commandfactory.hpp"

class GetCurrentCommand : public Command {

	public:
		GetCurrentCommand(string);
		GetCurrentCommand();
		~GetCurrentCommand() {}
	
		vector<string> setParameters();
		string getCommandName()			{ return "get.current";	}
		string getCommandCategory()		{ return "General";		}
        string getHelpString();	
        string getOutputPattern(string)	{ return ""; }
		string getCitation() { return "http://www.mothur.org/wiki/Get.current"; }
		string getDescription()		{ return "get current files saved by mothur"; }

	
		int execute(); 
		void help() { m->mothurOut(getHelpString()); }	
	
	
	private:
		
        CommandFactory* cFactory;
		vector<string> outputNames;
		bool abort;
	
		string clearTypes;
		vector<string> types;
		
};

#endif

