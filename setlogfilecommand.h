#ifndef SETLOGFILECOMMAND_H
#define SETLOGFILECOMMAND_H

/*
 *  setlogfilecommand.h
 *  Mothur
 *
 *  Created by westcott on 4/27/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "command.hpp"
#include "commandfactory.hpp"

/**********************************************************/

class SetLogFileCommand : public Command {
	
public:
	SetLogFileCommand(string);
	SetLogFileCommand() { setParameters(); abort = true; calledHelp = true; }
	~SetLogFileCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "set.logfile";		}
	string getCommandCategory()		{ return "General";			}
	string getHelpString();	
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	CommandFactory* commandFactory;
	string name;
	bool abort, append;
	vector<string> outputNames;
		
};

/**********************************************************/
 
#endif


