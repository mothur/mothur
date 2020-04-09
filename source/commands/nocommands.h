#ifndef NOCOMMAND_H
#define NOCOMMAND_H
/*
 *  nocommand.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

/* This command is run if the user enters an invalid command. */

#include "command.hpp"
#include "commandfactory.hpp"

class NoCommand : public Command {

public:
	NoCommand(string);
	~NoCommand(){}
	
	vector<string> setParameters()      { return outputNames;       } //dummy, doesn't really do anything	
	string getCommandName()             { return "NoCommand";       }
	string getCommandCategory()         { return "Hidden";          }
	string getHelpString()              { return "No Command";      }	
    string getOutputPattern(string) { return "";                }
	string getCitation()                { return "no citation";     }
	string getDescription()             { return "no description";  }

	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
	
private:
	vector<string> outputNames;
		
};

#endif
