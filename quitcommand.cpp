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
		
		else if (option != "") { cout << "There are no valid parameters for the quit command." << endl;  abort = true;  }

}
//**********************************************************************************************************************

void QuitCommand::help(){
	try {
		cout << "The quit command will terminate mothur and should be in the following format: " << "\n";
		cout << "quit() or quit" << "\n" << "\n";
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the QuitCommand class Function help. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the QuitCommand class function help. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
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
