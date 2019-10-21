#ifndef ENGINE_HPP
#define ENGINE_HPP

/*
 *  engine.hpp
 *  
 *
 *  Created by Pat Schloss on 8/15/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 */
 


#include "mothur.h"
#include "commandoptionparser.hpp"
#include "command.hpp"
#include "commandfactory.hpp"
#include "mothurout.h"

class Engine {
public:
	Engine(string);
	virtual ~Engine(){}
	virtual bool getInput() = 0;
	virtual string getLogFileName()			{	return mout->getLogFileName();          }

	vector<string> getOptions()		{	return options;		}
protected:
	vector<string> options;
	CommandFactory* cFactory;
	MothurOut* mout;
    CurrentFile* current;
    Utils util;
    time_t start;
    int numCommandsRun;
    string path;
};



class BatchEngine : public Engine {
public:
	BatchEngine(string, string);
	~BatchEngine();
	virtual bool getInput();
	bool openedBatch;
private:
	ifstream inputBatchFile;
	string getNextCommand(ifstream&);
    string batchFileName;

};



class InteractEngine : public Engine {
public:
	InteractEngine(string);
	~InteractEngine();
	virtual bool getInput();
    
private:
    string getCommand();
	
};


class ScriptEngine : public Engine {
public:
	ScriptEngine(string, string);
	~ScriptEngine();
	virtual bool getInput();
	bool openedBatch;
private:
	string listOfCommands;
	string getNextCommand(string&);

};


#endif
