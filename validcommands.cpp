/*
 *  validcommands.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/5/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "validcommands.h"

/***********************************************************************/

ValidCommands::ValidCommands() {
	try {
	
		commands["read.dist"]			= "read.dist"; 
		commands["read.otu"]			= "read.otu"; 
		commands["read.list"]			= "read.list"; 
		commands["cluster"]				= "cluster"; 
		commands["help"]				= "help"; 
		commands["quit"]				= "quit"; 
		commands["collect.single"]		= "collect.single"; 
		commands["collect.shared"]		= "collect.shared"; 
		commands["rarefaction.single"]	= "rarefaction.single"; 
		commands["rarefaction.shared"]	= "rarefaction.shared"; 
		commands["summary.single"]		= "summary.single"; 
		commands["summary.shared"]		= "summary.shared"; 
				
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ValidCommands class Function ValidCommands. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ValidCommands class function ValidCommands. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

/***********************************************************************/

ValidCommands::~ValidCommands() {}

/***********************************************************************/
bool ValidCommands::isValidCommand(string command) {
	try {	
	
		//is the command in the map
		if ((commands.find(command)) != (commands.end())) {
			return true;
		}else{
			cout << command << " is not a valid command in Mothur.  Valid commands are read.dist(), read.list(), cluster(), collect.single(), collect.shared(), rarefaction.single(), rarefaction.shared(), summary.single(), summary.shared(), quit(), help()." << endl;
			return false;
		}
		
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ValidCommands class Function isValidCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ValidCommands class function isValidCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/***********************************************************************/