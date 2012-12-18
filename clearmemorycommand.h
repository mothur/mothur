#ifndef CLEARMEMORYCOMMAND_H
#define CLEARMEMORYCOMMAND_H

/*
 *  clearmemorycommand.h
 *  Mothur
 *
 *  Created by westcott on 7/6/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "command.hpp"

class ClearMemoryCommand : public Command {
public:
	ClearMemoryCommand(string);
	ClearMemoryCommand(){ abort = true; calledHelp = true; }
	~ClearMemoryCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "clear.memory";			}
	string getCommandCategory()		{ return "General";	}
	
	string getHelpString();	
    string getOutputPattern(string) { return ""; }	
	string getCitation() { return "http://www.mothur.org/wiki/Clear.memory"; }
	string getDescription()		{ return "remove saved references from memory"; }
	
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
	
private:
	bool abort;
	vector<string> outputNames;
};

#endif

