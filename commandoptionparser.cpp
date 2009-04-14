/*
 *  commandoptionparser.cpp
 *  
 *
 *  Created by Pat Schloss on 10/23/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 */


using namespace std;


#include "globaldata.hpp"
#include "commandoptionparser.hpp"


//**********************************************************************************************************************
//This Function parses through the command line and pulls out the command then sends the options to  the parseGlobalData
CommandOptionParser::CommandOptionParser(string input){
	try {
		int openParen = input.find_first_of('(');
		int closeParen = input.find_last_of(')');
		string optionString = "";
		commandString = "";
	
		if(openParen != -1 && closeParen != -1){			
			commandString = input.substr(0, openParen);   //commandString contains everything before "("
			optionString = input.substr(openParen+1, closeParen-openParen-1); //optionString contains everything between "(" and ")".
		}
					
		GlobalData* globaldata = GlobalData::getInstance();
		globaldata->parseGlobalData(commandString, optionString);			//parser to separate and check options
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the CommandOptionParser class Function CommandOptionParser. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the CommandOptionParser class function CommandOptionParser. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}

}

//**********************************************************************************************************************

string CommandOptionParser::getCommandString()	{	return commandString;	}

//**********************************************************************************************************************
