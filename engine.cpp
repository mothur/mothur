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
	string logFileName = "mothur.logFile";
	remove(logFileName.c_str());
	
	system("clear");
}

/***********************************************************************/

InteractEngine::~InteractEngine(){
	}

/***********************************************************************/
//This function allows the user to input commands one line at a time until they quit.
//If the command is garbage it does nothing.
bool InteractEngine::getInput(){
	try {
		string input = "";
		string commandName = "";
		string options = "";
		int quitCommandCalled = 0;
		
		mothurOut("mothur v.1.4.1");
		mothurOutEndLine();		
		mothurOut("Last updated: 6/23/2009");
		mothurOutEndLine();	
		mothurOutEndLine();		
		mothurOut("by");
		mothurOutEndLine();		
		mothurOut("Patrick D. Schloss");
		mothurOutEndLine();
		mothurOutEndLine();			
		mothurOut("Department of Microbiology");
		mothurOutEndLine();		
		mothurOut("pschloss@micro.umass.edu");
		mothurOutEndLine();		
		mothurOut("http://schloss.micro.umass.edu/mothur");
		mothurOutEndLine();	
		mothurOutEndLine();	
		mothurOutEndLine();		
		mothurOut("Distributed under the GNU General Public License");
		mothurOutEndLine();
		mothurOutEndLine();			
		mothurOut("Type 'help()' for information on the commands that are available");
		mothurOutEndLine();
		mothurOutEndLine();			
		mothurOut("Type 'quit()' to exit program");
		mothurOutEndLine();	
		
		while(quitCommandCalled != 1){
			
			mothurOutEndLine();
			mothurOut("mothur > ");
			getline(cin, input);
			if (cin.eof()) { input = "quit()"; }
			
			mothurOutJustToLog(input);
			mothurOutEndLine();
			
			//allow user to omit the () on the quit command
			if (input == "quit") { input = "quit()"; }
			
			CommandOptionParser parser(input);
			commandName = parser.getCommandString();
			options = parser.getOptionString();
			
			if (commandName != "") {
			
				//executes valid command
				CommandFactory cFactory;
				Command* command = cFactory.getCommand(commandName, options);
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
//This function opens the batchfile to be used by BatchEngine::getInput.
BatchEngine::BatchEngine(string path, string batchFileName){
	try {
		globaldata = GlobalData::getInstance();
	
		openedBatch = openInputFile(batchFileName, inputBatchFile);
		globaldata->argv = path;
		string logFileName = "mothur.logFile";
		remove(logFileName.c_str());
		
		system("clear");
	
	//	char buffer = ' ';
	//	ifstream header("introtext.txt");
	//	while(!header.eof()){
	//		cout << buffer;
	//		buffer = header.get();
	//	}
	}
	catch(exception& e) {
		errorOut(e, "BatchEngine", "BatchEngine");
		exit(1);
	}
}

/***********************************************************************/

BatchEngine::~BatchEngine(){
	}

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
			else { getline(inputBatchFile, input); }
			
			
			
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
					CommandFactory cFactory;
					Command* command = cFactory.getCommand(commandName, options);
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
		string logFileName = "mothur.logFile";
		remove(logFileName.c_str());

		system("clear");
	
	}
	catch(exception& e) {
		errorOut(e, "ScriptEngine", "ScriptEngine");
		exit(1);
	}
}

/***********************************************************************/

ScriptEngine::~ScriptEngine(){
	}

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
				CommandFactory cFactory;
				Command* command = cFactory.getCommand(commandName, options);
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
		
		nextcommand = commandString.substr(0,commandString.find_first_of(';'));

				
		if ((commandString.find_first_of(';')+1) <= commandString.length()) {
			commandString = commandString.substr(commandString.find_first_of(';')+1, commandString.length());
		}else { commandString = ""; } //you have reached the last command.
		
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
