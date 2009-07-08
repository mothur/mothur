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
			//valid paramters for this command
			string Array[] =  {"command"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			
			//check to make sure all parameters are valid for command
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//check for required parameters
			command = validParameter.validFile(parameters, "command", false);
			if (command == "not found") { mothurOut("command is a required parameter."); mothurOutEndLine(); abort = true; }
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
		mothurOut("The system command parameter is command and it is required.\n");
		mothurOut("The system command should be in the following format: system(command=yourCommand).\n");
		mothurOut("Example system(command=clear).\n");
		mothurOut("Note: No spaces between parameter labels (i.e. command), '=' and parameters (i.e.yourCommand).\n\n");
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
