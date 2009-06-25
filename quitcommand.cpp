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

QuitCommand::QuitCommand(string option){
		abort = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else if (option != "") { mothurOut("There are no valid parameters for the quit command."); mothurOutEndLine();  abort = true;  }

}
//**********************************************************************************************************************

void QuitCommand::help(){
	try {
		 mothurOut("The quit command will terminate mothur and should be in the following format: \n"); 
		 mothurOut("quit() or quit\n\n");
	}
	catch(exception& e) {
		errorOut(e, "QuitCommand", "help");
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
