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

#include "mothurout.h"
#include "currentfile.h"

class Command;

class CommandFactory {
public:
	static CommandFactory* getInstance();
	Command* getCommand(string, string, string);
	Command* getCommand(string, string);
	Command* getCommand(string);

	bool isValidCommand(string);
	bool isValidCommand(string, string);
	void printCommands(ostream&);
    void printCommandsCategories(ostream&);
	map<string, string> getListCommands()	{	return commands;		}
	
private:
	MothurOut* m;
	CurrentFile* current;
    Utils util;
	
	map<string, string> commands;
	map<string, string>::iterator it;
	bool append;
	
    int checkForRedirects(string);
    
	static CommandFactory* _uniqueInstance;
	CommandFactory( const CommandFactory& ); // Disable copy constructor
	void operator=( const CommandFactory& ); // Disable assignment operator
	CommandFactory();
	~CommandFactory();

};

#endif
