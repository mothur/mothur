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
vector<string> NoCommand::getValidParameters(){	
	try {
		vector<string> myArray; 
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "NoCommand", "getValidParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> NoCommand::getRequiredParameters(){	
	try {
		vector<string> myArray;
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "NoCommand", "getRequiredParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> NoCommand::getRequiredFiles(){	
	try {
		vector<string> myArray;
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "NoCommand", "getRequiredFiles");
		exit(1);
	}
}
//**********************************************************************************************************************

NoCommand::NoCommand(string option)  {}

//**********************************************************************************************************************

NoCommand::~NoCommand(){}

//**********************************************************************************************************************

int NoCommand::execute(){
	//Could choose to give more help here?fdsah
	cout << "Invalid command.\n";
   
		CommandFactory* valid =  CommandFactory::getInstance();
	valid->printCommands(cout);
	
	return 0;
}

//**********************************************************************************************************************
