/*
 *  nocommand.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "nocommands.h"

//**********************************************************************************************************************

NoCommand::NoCommand(string option) : Command()  {}

//**********************************************************************************************************************

int NoCommand::execute(){
    MothurOut* m = MothurOut::getInstance();
	//Could choose to give more help here?
	m->mothurOut("[ERROR]: Invalid command.\n");
   
	CommandFactory* valid =  CommandFactory::getInstance();
	valid->printCommands(cout);
	
	return 2;
}

//**********************************************************************************************************************
