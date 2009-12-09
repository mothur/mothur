#ifndef COMMANDFACTORY_HPP
#define COMMANDFACTORY_HPP

/*
 *  commandfactory.h
 *  
 *
 *  Created by Pat Schloss on 10/25/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 */

#include "mothur.h"

class Command;

class CommandFactory {
public:
	static CommandFactory* getInstance();
	Command* getCommand(string, string);
	Command* getCommand();
	bool isValidCommand(string);
	void printCommands(ostream&);

private:
	Command* command;
	map<string, string> commands;
	map<string, string>::iterator it;
	
	static CommandFactory* _uniqueInstance;
	CommandFactory( const CommandFactory& ); // Disable copy constructor
	void operator=( const CommandFactory& ); // Disable assignment operator
	CommandFactory();
	~CommandFactory();

};

#endif
