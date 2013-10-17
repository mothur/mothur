/*
 *  commandoptionparser.cpp
 *  
 *
 *  Created by Pat Schloss on 10/23/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 */


#include "commandoptionparser.hpp"


//**********************************************************************************************************************
//This Function parses through the command line and pulls out the command then sends the options to  the parseGlobalData
CommandOptionParser::CommandOptionParser(string input){
	try {
		m = MothurOut::getInstance();
		
		int openParen = input.find_first_of('(');
		int closeParen = input.find_last_of(')');
		optionString = "";
		commandString = "";

		if(openParen != -1 && closeParen != -1){	
            //gobble extra spaces
            int spot = 0;
            for (int i = 0; i < input.length(); i++) {  if (!(isspace(input[i]))) { spot = i; break; } }
            if (spot > openParen) { spot = 0; }
			commandString = input.substr(spot, openParen-spot);   //commandString contains everything before "("
			optionString = input.substr((openParen+1), (closeParen-openParen-1)); //optionString contains everything between "(" and ")".
		}
		else if (openParen == -1) { m->mothurOut("[ERROR]: You are missing ("); m->mothurOutEndLine(); }
		else if (closeParen == -1) { m->mothurOut("[ERROR]: You are missing )"); m->mothurOutEndLine(); }
    }
	catch(exception& e) {
		m->errorOut(e, "CommandOptionParser", "CommandOptionParser");
		exit(1);
	}
}

//**********************************************************************************************************************

string CommandOptionParser::getCommandString()	{	return commandString;	}

//**********************************************************************************************************************

string CommandOptionParser::getOptionString()	{	return optionString;	}

//**********************************************************************************************************************
