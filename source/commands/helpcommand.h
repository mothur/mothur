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
	~HelpCommand(){}
	
	vector<string> setParameters()	{ return outputNames;	} //dummy, doesn't really do anything	
	string getCommandName()			{ return "help";		}
	string getCommandCategory()		{ return "Hidden";		}
	string getHelpString() { return "For more information about a specific command type 'commandName(help)' i.e. 'cluster(help)'"; }
    string getCommonQuestions();
    string getOutputPattern(string) { return "";                }
	string getCitation() { return "no citation"; }
	string getDescription()		{ return "help"; }

	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
	
private:
	CommandFactory* validCommands;
	vector<string> outputNames;
    
    bool abort, calledHelp;
    
    string commandName;
		
};
 
#endif
