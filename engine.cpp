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
Engine::Engine(){
	try {
		cFactory = CommandFactory::getInstance();
		mout = MothurOut::getInstance();
	}
	catch(exception& e) {
		mout->errorOut(e, "Engine", "Engine");
		exit(1);
	}
}

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

			mout->mothurOutEndLine();
			
			input = getCommand();	
			mout->mothurOutEndLine();	
			
			if (mout->control_pressed) { input = "quit()"; }
			
			//allow user to omit the () on the quit command
			if (input == "quit") { input = "quit()"; }
			
			CommandOptionParser parser(input);
			commandName = parser.getCommandString();
	
			options = parser.getOptionString();
			
			if (commandName != "") {
				mout->executing = true;
				
				#ifdef USE_MPI
					int pid;
					MPI_Comm_rank(MPI_COMM_WORLD, &pid); 
					
					MPI_Barrier(MPI_COMM_WORLD); //make everyone wait - just in case
				
					if ((cFactory->MPIEnabled(commandName)) || (pid == 0)) {
				#endif
				//executes valid command
				Command* command = cFactory->getCommand(commandName, options);
				quitCommandCalled = command->execute();
				mout->control_pressed = 0;
				mout->executing = false;
				
				#ifdef USE_MPI
					}
				#endif
			}else {
				mout->mothurOut("Your input contains errors. Please try again."); 
				mout->mothurOutEndLine();
			}
		}	
		return 1;
	}
	catch(exception& e) {
		mout->errorOut(e, "InteractEngine", "getInput");
		exit(1);
	}
}
/***********************************************************************/
string Engine::getCommand()  {
	try {
		#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
			#ifdef USE_READLINE
				char* nextCommand = NULL;
				nextCommand = readline("mothur > ");
				
				if(nextCommand != NULL) {  add_history(nextCommand);  }	
				else{ //^D causes null string and we want it to quit mothur
					nextCommand = "quit"; 
					mout->mothurOut(nextCommand);
				}	
				
				mout->mothurOutJustToLog("mothur > " + toString(nextCommand));
				return nextCommand;
			#else
				string nextCommand = "";
				mout->mothurOut("mothur > ");
				getline(cin, nextCommand);
				mout->mothurOutJustToLog("mothur > " + toString(nextCommand));
				
				return nextCommand;
			#endif
		#else
				string nextCommand = "";
				
				mout->mothurOut("mothur > ");
				getline(cin, nextCommand);
				mout->mothurOutJustToLog(toString(nextCommand));
				
				return nextCommand;
		#endif
	
						
	}
	catch(exception& e) {
		mout->errorOut(e, "Engine", "getCommand");
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
		mout->errorOut(e, "BatchEngine", "BatchEngine");
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
			mout->mothurOut("unable to open batchfile");  
			mout->mothurOutEndLine();
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
				
				mout->mothurOutEndLine();
				mout->mothurOut("mothur > " + input);
				mout->mothurOutEndLine();
							
				if (mout->control_pressed) { input = "quit()"; }
				
				//allow user to omit the () on the quit command
				if (input == "quit") { input = "quit()"; }

				CommandOptionParser parser(input);
				commandName = parser.getCommandString();
				options = parser.getOptionString();
										
				if (commandName != "") {
					mout->executing = true;
					#ifdef USE_MPI
						int pid;
						MPI_Comm_rank(MPI_COMM_WORLD, &pid); 
cout << pid << " is waiting " << commandName << endl;						
						MPI_Barrier(MPI_COMM_WORLD); //make everyone wait - just in case
cout << pid << " is here " << commandName << endl;
						if ((cFactory->MPIEnabled(commandName)) || (pid == 0)) {
					#endif
					//executes valid command
					Command* command = cFactory->getCommand(commandName, options);
					quitCommandCalled = command->execute();
					mout->control_pressed = 0;
					mout->executing = false;
				
					#ifdef USE_MPI
						}
					#endif
				}else {		
					mout->mothurOut("Invalid."); 
					mout->mothurOutEndLine();
				}
				
			}
			gobble(inputBatchFile);
		}
		
		inputBatchFile.close();
		return 1;
	}
	catch(exception& e) {
		mout->errorOut(e, "BatchEngine", "getInput");
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
		mout->errorOut(e, "ScriptEngine", "ScriptEngine");
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
			
			mout->mothurOutEndLine();
			mout->mothurOut("mothur > " + input);
			mout->mothurOutEndLine();
			
			if (mout->control_pressed) { input = "quit()"; }
				
			//allow user to omit the () on the quit command
			if (input == "quit") { input = "quit()"; }

			CommandOptionParser parser(input);
			commandName = parser.getCommandString();
			options = parser.getOptionString();
										
			if (commandName != "") {
				mout->executing = true;
				#ifdef USE_MPI
					int pid;
					MPI_Comm_rank(MPI_COMM_WORLD, &pid); 
					
					MPI_Barrier(MPI_COMM_WORLD); //make everyone wait - just in case
					
					if ((cFactory->MPIEnabled(commandName)) || (pid == 0)) {
				#endif
				//executes valid command
				Command* command = cFactory->getCommand(commandName, options);
				quitCommandCalled = command->execute();
				mout->control_pressed = 0;
				mout->executing = false;
				
				#ifdef USE_MPI
					}
				#endif
			}else {		
				mout->mothurOut("Invalid."); 
				mout->mothurOutEndLine();
			}
			
		}
		
		return 1;
	}
	catch(exception& e) {
		mout->errorOut(e, "ScriptEngine", "getInput");
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
		mout->errorOut(e, "ScriptEngine", "getNextCommand");
		exit(1);
	}
}
/***********************************************************************/
