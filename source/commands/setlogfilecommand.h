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
	~SetLogFileCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "set.logfile";		}
	string getCommandCategory()		{ return "General";			}
	
	string getHelpString();	
    string getOutputPattern(string){ return ""; }	
	string getCitation() { return "http://www.mothur.org/wiki/Set.logfile"; }
	string getDescription()		{ return "set logfile name"; }

	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	string name;
	bool abort, append;
	vector<string> outputNames;
		
};

/**********************************************************/
 
#endif


