/*
 *  engine.cpp
 *  
 *
 *  Created by Pat Schloss on 8/15/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 *  There's a TON of duplicated code between InteractEngine and BatchEngine
 *  I couldn't figure out how to transition between ifstream (batch) and cin (interact)
 *  Fix later, don't have time now.
 *
 */


#include "engine.hpp"


/***********************************************************************/

InteractEngine::InteractEngine(string path){

	globaldata = GlobalData::getInstance();
	globaldata->argv = path;
}

/***********************************************************************/

InteractEngine::~InteractEngine(){}

/***********************************************************************/
//This function allows the user to input commands one line at a time until they quit.
//If the command is garbage it does nothing.
bool InteractEngine::getInput(){
	try {
		string input = "";
		string commandName = "";
		string options = "";
		int quitCommandCalled = 0;
		
		while(quitCommandCalled != 1){
			
			mothurOutEndLine();
			
			input = getCommand();			
			
			mothurOutJustToLog("mothur > " + input);
			mothurOutEndLine();
			
			//allow user to omit the () on the quit command
			if (input == "quit") { input = "quit()"; }
			
			CommandOptionParser parser(input);
			commandName = parser.getCommandString();
	
			options = parser.getOptionString();
			
			if (commandName != "") {
			
				//executes valid command
				Command* command = cFactory->getCommand(commandName, options);
				quitCommandCalled = command->execute();
				
			}else {
				mothurOut("Your input contains errors. Please try again."); 
				mothurOutEndLine();
			}
		}	
		return 1;
	}
	catch(exception& e) {
		errorOut(e, "InteractEngine", "getInput");
		exit(1);
	}
}
/***********************************************************************/
string Engine::getCommand()  {
	try {
		char* nextCommand = NULL;
		
		nextCommand = readline("mothur > ");
		if(nextCommand != NULL) {  add_history(nextCommand);  }
						
		return nextCommand;
	}
	catch(exception& e) {
		errorOut(e, "Engine", "getCommand");
		exit(1);
	}
}
/***********************************************************************/
void Engine::terminateCommand(int dummy)  {
	try {
		mothurOut("Stopping command."); mothurOutEndLine();
		cFactory->getCommand();  //terminates command
	}
	catch(exception& e) {
		errorOut(e, "Engine", "terminateCommand");
		exit(1);
	}
}
/***********************************************************************/
//This function opens the batchfile to be used by BatchEngine::getInput.
BatchEngine::BatchEngine(string path, string batchFileName){
	try {
		globaldata = GlobalData::getInstance();
	
		openedBatch = openInputFile(batchFileName, inputBatchFile);
		globaldata->argv = path;
				
	}
	catch(exception& e) {
		errorOut(e, "BatchEngine", "BatchEngine");
		exit(1);
	}
}

/***********************************************************************/

BatchEngine::~BatchEngine(){	}

/***********************************************************************/
//This Function allows the user to run a batchfile containing several commands on Dotur
bool BatchEngine::getInput(){
	try {
		//check if this is a valid batchfile
		if (openedBatch == 1) {  
			mothurOut("unable to open batchfile");  
			mothurOutEndLine();
			return 1; 
		}
	
		string input = "";
		string commandName = "";
		string options = "";
		
		//CommandFactory cFactory;
		int quitCommandCalled = 0;
	
		while(quitCommandCalled == 0){
	
			if (inputBatchFile.eof()) { input = "quit()"; }
			else { input = getline(inputBatchFile); }
			
			if (input[0] != '#') {
				
				mothurOutEndLine();
				mothurOut("mothur > " + input);
				mothurOutEndLine();
				
				
				//allow user to omit the () on the quit command
				if (input == "quit") { input = "quit()"; }

				CommandOptionParser parser(input);
				commandName = parser.getCommandString();
				options = parser.getOptionString();
										
				if (commandName != "") {

					//executes valid command
					Command* command = cFactory->getCommand(commandName, options);
					quitCommandCalled = command->execute();
				}else {		
					mothurOut("Invalid."); 
					mothurOutEndLine();
				}
				
			}
			gobble(inputBatchFile);
		}
		
		inputBatchFile.close();
		return 1;
	}
	catch(exception& e) {
		errorOut(e, "BatchEngine", "getInput");
		exit(1);
	}
}


/***********************************************************************/
/***********************************************************************/
//This function opens the batchfile to be used by BatchEngine::getInput.
ScriptEngine::ScriptEngine(string path, string commandString){
	try {
		globaldata = GlobalData::getInstance();
		
		//remove quotes
		listOfCommands = commandString.substr(1, (commandString.length()-1));

		globaldata->argv = path;
		
	}
	catch(exception& e) {
		errorOut(e, "ScriptEngine", "ScriptEngine");
		exit(1);
	}
}

/***********************************************************************/

ScriptEngine::~ScriptEngine(){ 	}

/***********************************************************************/
//This Function allows the user to run a batchfile containing several commands on mothur
bool ScriptEngine::getInput(){
	try {
			
		string input = "";
		string commandName = "";
		string options = "";
		
		
		//CommandFactory cFactory;
		int quitCommandCalled = 0;
	
		while(quitCommandCalled == 0){
		
			input = getNextCommand(listOfCommands);	
			
			if (input == "") { input = "quit()"; }
			
			
			mothurOutEndLine();
			mothurOut("mothur > " + input);
			mothurOutEndLine();

				
			//allow user to omit the () on the quit command
			if (input == "quit") { input = "quit()"; }

			CommandOptionParser parser(input);
			commandName = parser.getCommandString();
			options = parser.getOptionString();
										
			if (commandName != "") {

				//executes valid command
				Command* command = cFactory->getCommand(commandName, options);
				quitCommandCalled = command->execute();
			}else {		
				mothurOut("Invalid."); 
				mothurOutEndLine();
			}
			
		}
		
		return 1;
	}
	catch(exception& e) {
		errorOut(e, "ScriptEngine", "getInput");
		exit(1);
	}
}
/***********************************************************************/
string ScriptEngine::getNextCommand(string& commandString) {
	try {
		string nextcommand = "";
		int count = 0;
		
		//go through string until you reach ; or end
		while (count < commandString.length()) { 
			
			if (commandString[count] == ';') {  break;   }
			else {		nextcommand += commandString[count];	}
			
			count++;
		}
		
		//if you are not at the end
		if (count != commandString.length())  {   commandString = commandString.substr(count+1, commandString.length());  }
		else { commandString = ""; }
				
		
		//get rid of spaces in between commands if any
		if (commandString.length() > 0) {
			while (commandString[0] == ' ') {  
				commandString = commandString.substr(1,commandString.length());
				if (commandString.length() == 0) {  break;  }
			}
		}
					
		return nextcommand;
	}
	catch(exception& e) {
		errorOut(e, "ScriptEngine", "getNextCommand");
		exit(1);
	}
}
/***********************************************************************/
