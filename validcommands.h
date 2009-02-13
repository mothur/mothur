#ifndef VALIDCOMMANDS_H
#define VALIDCOMMANDS_H

/*
 *  validcommands.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/5/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */
using namespace std;

#include <string>
#include <iostream>
#include <map>

//This class contains a list of all valid commands in Mothur.  
//It has a function which will tell you if your command is valid.
//When adding a new command you must add it to the valid list in the class constructor.

class ValidCommands {

	public:
		ValidCommands();
		~ValidCommands();
		bool isValidCommand(string);
		void printCommands(ostream&);
	private:
		map<string, string> commands;
		map<string, string>::iterator it;

};

#endif
