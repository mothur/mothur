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

NoCommand::NoCommand(string option)  {}

//**********************************************************************************************************************

int NoCommand::execute(){
	//Could choose to give more help here?fdsah
	cout << "Invalid command.\n";
   
	CommandFactory* valid =  CommandFactory::getInstance();
	valid->printCommands(cout);
	
	return 0;
}

//**********************************************************************************************************************
