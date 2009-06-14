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
#include "validcommands.h"


class HelpCommand : public Command {
	
public:
	HelpCommand(string);
	~HelpCommand();
	int execute();
	void help() {};
	
private:
	ValidCommands* validCommands;
	
private:
		
};
 
#endif
