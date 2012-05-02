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
QuitCommand::QuitCommand(string option) {
		abort = false; calledHelp = false;   
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}

}
//**********************************************************************************************************************
QuitCommand::~QuitCommand(){}
//**********************************************************************************************************************
int QuitCommand::execute(){
	if (abort == true) { return 0; }
	return 1;
}
//**********************************************************************************************************************
