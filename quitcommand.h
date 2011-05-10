#ifndef QUITCOMMAND_H
#define QUITCOMMAND_H
/*
 *  quitcommand.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "command.hpp"

/* The quit() command:
	The quit command terminates the mothur program. 
	The quit command should be in the following format: quit ().   */


class QuitCommand : public Command {
	
public:
	QuitCommand(string);
	QuitCommand() {}
	~QuitCommand();
	
	vector<string> setParameters()	{ return outputNames;	} //dummy, doesn't really do anything	
	string getCommandName()			{ return "quit";		}
	string getCommandCategory()		{ return "Hidden";		}
	string getHelpString() { return "The quit command will terminate mothur and should be in the following format: quit() or quit. \n"; }	
	string getCitation() { return "no citation"; }
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
	
private:
	bool abort;
	vector<string> outputNames;
};

#endif
