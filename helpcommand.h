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

#include <Carbon/Carbon.h>
#include "command.hpp"
#include "globaldata.hpp"


class HelpCommand : public Command {
	
public:
	HelpCommand();
	~HelpCommand();
	int execute();
private:
	GlobalData* globaldata;
	
private:
		
};
 
#endif