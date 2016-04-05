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
#include "mothurout.h"
#include "currentfile.h"

class Command;

class CommandFactory {
public:
	static CommandFactory* getInstance();
	Command* getCommand(string, string, string);
	Command* getCommand(string, string);
	Command* getCommand(string);
	//Command* getCommand();
	bool isValidCommand(string);
	bool isValidCommand(string, string);
	void printCommands(ostream&);
    void printCommandsCategories(ostream&);
	void setOutputDirectory(string o)		{	if(m->dirCheck(o) || (o == "")) {  outputDir = o; m->setOutputDir(o); }	}
	void setInputDirectory(string i)		{	if(m->dirCheck(i) || (i == "")) {  inputDir = i;	}	}
	void setLogfileName(string n, bool a)	{	logFileName = n;  append = a;		}
	string getLogfileName()					{	return logFileName; 	}
	bool getAppend()						{	return append;			}
	string getOutputDir()					{	return outputDir;		}
    string getInputDir()					{	return inputDir;		}
	map<string, string> getListCommands()	{	return commands;		}
	
private:
	Command* command;
	Command* shellcommand;
	Command* pipecommand;
	
	MothurOut* m;
	CurrentFile* currentFile;
	
	map<string, string> commands;
	map<string, string>::iterator it;
	string outputDir, inputDir, logFileName;
	bool append;
	
    int checkForRedirects(string);
    
	static CommandFactory* _uniqueInstance;
	CommandFactory( const CommandFactory& ); // Disable copy constructor
	void operator=( const CommandFactory& ); // Disable assignment operator
	CommandFactory();
	~CommandFactory();

};

#endif
