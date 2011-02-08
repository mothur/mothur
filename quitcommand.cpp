/*
 *  quitcommand.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "quitcommand.h"

//**********************************************************************************************************************
vector<string> QuitCommand::getValidParameters(){	
	try {
		vector<string> myArray; 
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "QuitCommand", "getValidParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> QuitCommand::getRequiredParameters(){	
	try {
		vector<string> myArray;
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "QuitCommand", "getRequiredParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> QuitCommand::getRequiredFiles(){	
	try {
		vector<string> myArray;
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "QuitCommand", "getRequiredFiles");
		exit(1);
	}
}
//**********************************************************************************************************************

QuitCommand::QuitCommand(string option) {
		abort = false; calledHelp = false;   
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }

}
//**********************************************************************************************************************

void QuitCommand::help(){
	try {
		 m->mothurOut("The quit command will terminate mothur and should be in the following format: \n"); 
		 m->mothurOut("quit() or quit\n\n");
	}
	catch(exception& e) {
		m->errorOut(e, "QuitCommand", "help");
		exit(1);
	}
}

//**********************************************************************************************************************

QuitCommand::~QuitCommand(){}

//**********************************************************************************************************************

int QuitCommand::execute(){
	if (abort == true) { return 0; }
	return 1;
}

//**********************************************************************************************************************
