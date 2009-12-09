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
#include "globaldata.hpp"
#include "commandoptionparser.hpp"
#include "command.hpp"
#include "commandfactory.hpp"

class GlobalData;

class Engine {
public:
	Engine(); 
	virtual ~Engine(){}
	virtual bool getInput() = 0;
	virtual string getCommand();
	vector<string> getOptions() {	return options;		}
	//virtual void terminateCommand(int);
protected:
	vector<string> options;
	CommandFactory* cFactory;
};



class BatchEngine : public Engine {
public:
	BatchEngine(string, string);
	~BatchEngine();
	virtual bool getInput();
	int openedBatch;
private:
	GlobalData* globaldata;
	ifstream inputBatchFile;

};



class InteractEngine : public Engine {
public:
	InteractEngine(string);
	~InteractEngine();
	virtual bool getInput();
private:
	GlobalData* globaldata;
};


class ScriptEngine : public Engine {
public:
	ScriptEngine(string, string);
	~ScriptEngine();
	virtual bool getInput();
	int openedBatch;
private:
	GlobalData* globaldata;
	string listOfCommands;
	string getNextCommand(string&);

};


#endif
