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

NoCommand::NoCommand(string option){}

//**********************************************************************************************************************

NoCommand::~NoCommand(){}

//**********************************************************************************************************************

int NoCommand::execute(){
	//Could choose to give more help here?fdsah
	cout << "Invalid command." << "\n";
	
	CommandFactory* valid = new CommandFactory();
	valid->printCommands(cout);
	delete valid;
	
	return 0;
}

//**********************************************************************************************************************
