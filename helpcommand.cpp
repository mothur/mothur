/*
 *  helpcommand.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "helpcommand.h"

//**********************************************************************************************************************
vector<string> HelpCommand::getValidParameters(){	
	try {
		
		vector<string> myArray; 
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "HelpCommand", "getValidParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> HelpCommand::getRequiredParameters(){	
	try {
		vector<string> myArray;
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "HelpCommand", "getRequiredParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> HelpCommand::getRequiredFiles(){	
	try {
		vector<string> myArray;
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "HelpCommand", "getRequiredFiles");
		exit(1);
	}
}
//**********************************************************************************************************************

HelpCommand::HelpCommand(string option)  {
	
	validCommands = CommandFactory::getInstance();
	
}

//**********************************************************************************************************************

HelpCommand::~HelpCommand(){}

//**********************************************************************************************************************

int HelpCommand::execute(){
	try {
		validCommands->printCommands(cout);
		m->mothurOut("For more information about a specific command type 'commandName(help)' i.e. 'read.dist(help)'"); m->mothurOutEndLine();
	
		m->mothurOutEndLine(); m->mothurOut("For further assistance please refer to the Mothur manual on our wiki at http://www.mothur.org/wiki, or contact Pat Schloss at mothur.bugs@gmail.com.\n");
	
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "HelpCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************/
