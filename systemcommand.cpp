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
SystemCommand::SystemCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		
		else {
			if (option == "") { m->mothurOut("You must enter a command to run."); m->mothurOutEndLine(); abort = true; }
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
		m->errorOut(e, "SystemCommand", "SystemCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

string SystemCommand::getHelpString(){
	try {
		string helpString = "";
		helpString += "The system command allows you to execute a system command from within mothur.\n";
		helpString += "The system has no parameters.\n";
		helpString += "The system command should be in the following format: system(yourCommand).\n";
		helpString += "Example system(clear).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "SystemCommand", "help");
		exit(1);
	}
}

//**********************************************************************************************************************

int SystemCommand::execute(){
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		system(command.c_str());
		
		return 0;		
	}

	catch(exception& e) {
		m->errorOut(e, "SystemCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************
