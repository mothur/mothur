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
vector<string> SystemCommand::setParameters(){	
	try {
		CommandParameter pcommand("command", "String", "", "", "", "", "",false,false); parameters.push_back(pcommand);
				
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "SystemCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
SystemCommand::SystemCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string, string> parameters = parser.getParameters();
			map<string, string>::iterator it;
			
			ValidParameters validParameter;
			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			string commandOption = validParameter.validFile(parameters, "command", false);			
			if (commandOption == "not found") { commandOption = ""; }
			else { command = commandOption; }
			
			if ((option == "") && (commandOption == "")) { m->mothurOut("You must enter a command to run."); m->mothurOutEndLine(); abort = true; }
			else if (commandOption == "") { 
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
