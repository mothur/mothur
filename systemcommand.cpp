/*
 *  systemcommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 7/8/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "systemcommand.h"

//**********************************************************************************************************************

SystemCommand::SystemCommand(string option){
	try {
		abort = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			if (option == "") { mothurOut("You must enter a command to run."); mothurOutEndLine(); abort = true; }
			else { 
				//check for outputdir and inputdir parameters
				int commaPos = option.find_first_of(',');
				
				//if there is a comma then grab string up to that pos
				if (commaPos != option.npos) {
					option = option.substr(0, commaPos);
				}
			
				command = option;
			}
		}	

	}
	catch(exception& e) {
		errorOut(e, "SystemCommand", "SystemCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

void SystemCommand::help(){
	try {
		mothurOut("The system command allows you to execute a system command from within mothur.\n");
		mothurOut("The system has no parameters.\n");
		mothurOut("The system command should be in the following format: system(yourCommand).\n");
		mothurOut("Example system(clear).\n");
	}
	catch(exception& e) {
		errorOut(e, "SystemCommand", "help");
		exit(1);
	}
}

//**********************************************************************************************************************

int SystemCommand::execute(){
	try {
		
		if (abort == true) { return 0; }
		
		system(command.c_str());
		
		return 0;		
	}

	catch(exception& e) {
		errorOut(e, "SystemCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************
