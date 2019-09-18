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
#include "mothur.h"

/***********************************************************************/
Engine::Engine(){
	try {
		cFactory = CommandFactory::getInstance();
		mout = MothurOut::getInstance();
        current = CurrentFile::getInstance();
        
        start = time(NULL);
        numCommandsRun = 0;
	}
	catch(exception& e) {
		mout->errorOut(e, "Engine", "Engine");
		exit(1);
	}
}
/***********************************************************************/

/***********************************************************************/

InteractEngine::InteractEngine(string path){

	
	string temppath = path.substr(0, (path.find_last_of("othur")-5));
	
	//this will happen if you set the path variable to contain mothur's exe location
    Utils util;
	if (temppath == "") { path = util.findProgramPath("mothur"); }
    else { path = temppath; }
	
	current->setProgramPath(util.getFullPathName(path));
    current->setBlastPath(current->getProgramPath());

    //if you haven't set your own location
    #ifdef MOTHUR_FILES
    #else
        //set default location to search for files to mothur's executable location.  This will resolve issue of double-clicking on the executable which opens mothur and sets pwd to your home directory instead of the mothur directory and leads to "unable to find file" errors.
        if (current->getProgramPath() != "") { current->setDefaultPath(current->getProgramPath()); }
    #endif
    
    if (mout->getLogFileName() == "") {
        time_t ltime = time(NULL); /* calendar time */
        string outputPath = current->getOutputDir();
        //if (outputPath == "") { outputPath = current->getDefaultPath();  }
        string logFileName = outputPath + "mothur." + toString(ltime) + ".logfile";
        mout->setLogFileName(logFileName, false);
        mout->mothurOut("\n");
    }
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
            
            input = getCommand();
            
            if (mout->getControl_pressed()) { input = "quit()"; }
            
            //allow user to omit the () on the quit command
            if (input == "quit") { input = "quit()"; }
            if (input == "help") { input = "help()"; }
            
            CommandOptionParser parser(input);
            commandName = parser.getCommandString();
            
            options = parser.getOptionString();
            
            if (commandName != "") {
                numCommandsRun++;
                mout->setExecuting(true);
                mout->resetCommandErrors();
                
                //executes valid command
                mout->setChangedSeqNames(true);
                
                Command* command = cFactory->getCommand(commandName, options);
                quitCommandCalled = command->execute();
                delete command;
                
                //if we aborted command
                if (quitCommandCalled == 2) {  mout->mothurOut("[ERROR]: did not complete " + commandName + ".\n");  }
                
                mout->setControl_pressed(false);
                mout->setExecuting(false);
                
            }else { mout->mothurOut("[ERROR]: Invalid.\n"); }
        }
		return 1;
	}
	catch(exception& e) {
		mout->errorOut(e, "InteractEngine", "getInput");
		exit(1);
	}
}
/***********************************************************************/
string InteractEngine::getCommand()  {
	try {
	
		#if defined NON_WINDOWS
			#ifdef USE_READLINE
				char* nextCommand = NULL;
				nextCommand = readline("\nmothur > ");
				
				if(nextCommand != NULL) {  add_history(nextCommand);  }	
				else{ //^D causes null string and we want it to quit mothur
					nextCommand = strdup("quit");
				}	
				
				mout->mothurOutJustToLog("\nmothur > " + toString(nextCommand) + "\n");
                string returnCommand = nextCommand;
                free(nextCommand);
				return returnCommand;
			#else
				string nextCommand = "";
				mout->mothurOut("\nmothur > ");
				getline(cin, nextCommand);
                mout->mothurOut("\n");
				mout->mothurOutJustToLog("\nmothur > " + toString(nextCommand) + "\n");
				
				return nextCommand;
			#endif
		#else
				string nextCommand = "";
				
				mout->mothurOut("\nmothur > ");
				getline(cin, nextCommand);
                mout->mothurOut("\n");
				mout->mothurOutJustToLog(toString(nextCommand) + "\n");
				
				return nextCommand;
		#endif
	
						
	}
	catch(exception& e) {
		mout->errorOut(e, "InteractEngine", "getCommand");
		exit(1);
	}
}
/***********************************************************************/
//This function opens the batchfile to be used by BatchEngine::getInput.
BatchEngine::BatchEngine(string path, string batchFileName){
	try {
		string temppath = path.substr(0, (path.find_last_of("othur")-5));
	
		//this will happen if you set the path variable to contain mothur's exe location
		if (temppath == "") { path = util.findProgramPath("mothur"); }
        else { path = temppath; }
		
        current->setProgramPath(util.getFullPathName(path));
        current->setBlastPath(current->getProgramPath());
        
        openedBatch = util.openInputFile(batchFileName, inputBatchFile, "no error");
        if (!openedBatch) {
            if (util.checkLocations(batchFileName, current->getLocations())) { openedBatch = util.openInputFile(batchFileName, inputBatchFile); }
            else {  mout->mothurOut("[ERROR]: unable to open batch file, please correct.\n");  }
        }
        
        //if you haven't set your own location
#ifdef MOTHUR_FILES
#else
        //set default location to search for files to mothur's executable location.  This will resolve issue of double-clicking on the executable which opens mothur and sets pwd to your home directory instead of the mothur directory and leads to "unable to find file" errors.
        if (current->getProgramPath() != "") { current->setDefaultPath(current->getProgramPath()); }
#endif

				
	}
	catch(exception& e) {
		mout->errorOut(e, "BatchEngine", "BatchEngine");
		exit(1);
	}
}

/***********************************************************************/

BatchEngine::~BatchEngine(){
    time_t end = time(NULL);
    mout->mothurOut("\n\nIt took " + toString(end-start) + " seconds to run " + toString(numCommandsRun) + " commands from your batch file.\n\n");
}

/***********************************************************************/
//This Function allows the user to run a batchfile containing several commands on Dotur
bool BatchEngine::getInput(){
	try {
		//check if this is a valid batchfile
		if (!openedBatch) { return 1;  }
	
		string input = "";
		string commandName = "";
		string options = "";
		
		int quitCommandCalled = 0;
	    int count = 0;
		while(quitCommandCalled != 1){
			
			input = getNextCommand(inputBatchFile);
			count++;
			
            if (input[0] != '#') {
				mout->appendLogBuffer("\nmothur > " + input + "\n");
				
				if (mout->getControl_pressed()) { input = "quit()"; }
				
				//allow user to omit the () on the quit command
				if (input == "quit") { input = "quit()"; }
                if (input == "help") { input = "help()"; }

				CommandOptionParser parser(input);
				commandName = parser.getCommandString();
				options = parser.getOptionString();
										
				if (commandName != "") {
                    numCommandsRun++;
					mout->setExecuting(true);
                    mout->resetCommandErrors();
					
					//executes valid command
                    mout->setChangedSeqNames(true);
							
					Command* command = cFactory->getCommand(commandName, options);
					quitCommandCalled = command->execute();
                    delete command;
							
					//if we aborted command
					if (quitCommandCalled == 2) {  mout->mothurOut("[ERROR]: did not complete " + commandName + ".\n");  }
                    
                    if (mout->getControl_pressed()) { break;  }
                    mout->setControl_pressed(false);
                    mout->setExecuting(false);
										
				}else {	 mout->mothurOut("[ERROR]: Invalid.\n"); }
				
			}
			util.gobble(inputBatchFile);
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
string BatchEngine::getNextCommand(ifstream& inputBatchFile) {
	try {
			
		string nextcommand = "";
		
		if (inputBatchFile.eof()) { nextcommand = "quit()"; }
		else { nextcommand = util.getline(inputBatchFile); }
		
		return nextcommand;
	}
	catch(exception& e) {
		mout->errorOut(e, "BatchEngine", "getNextCommand");
		exit(1);
	}
}

/***********************************************************************/
/***********************************************************************/
ScriptEngine::ScriptEngine(string path, string commandString){
	try {
		
		//remove quotes
		listOfCommands = commandString.substr(1, (commandString.length()-1));
		
		string temppath = path.substr(0, (path.find_last_of("othur")-5));

		//this will happen if you set the path variable to contain mothur's exe location
		if (temppath == "") { path = util.findProgramPath("mothur"); }
        else { path = temppath; }
		
        current->setProgramPath(util.getFullPathName(path));
        current->setBlastPath(current->getProgramPath());
    
        //if you haven't set your own location
#ifdef MOTHUR_FILES
#else
        //set default location to search for files to mothur's executable location.  This will resolve issue of double-clicking on the executable which opens mothur and sets pwd to your home directory instead of the mothur directory and leads to "unable to find file" errors.
        //string tempDefault = path.substr(0, (path.find_last_of('m')));
        if (path != "") { current->setDefaultPath(path); }
#endif

				
	}
	catch(exception& e) {
		mout->errorOut(e, "ScriptEngine", "ScriptEngine");
		exit(1);
	}
}

/***********************************************************************/

ScriptEngine::~ScriptEngine(){
    time_t end = time(NULL);
    mout->mothurOut("\n\nIt took " + toString(end-start) + " seconds to run " + toString(numCommandsRun) + " commands from your script.\n\n");
}

/***********************************************************************/
//This Function allows the user to run a batchfile containing several commands on mothur
bool ScriptEngine::getInput(){
	try {
			
		string input = "";
		string commandName = "";
		string options = "";
		
		int quitCommandCalled = 0;
	
		while(quitCommandCalled != 1){
			
			input = getNextCommand(listOfCommands);	
			
			if (input == "") { input = "quit()"; }
			
            mout->appendLogBuffer("\nmothur > " + input + "\n");
        
            if (mout->getControl_pressed()) { input = "quit()"; }
				
			//allow user to omit the () on the quit command
			if (input == "quit") { input = "quit()"; }
            if (input == "help") { input = "help()"; }

			CommandOptionParser parser(input);
			commandName = parser.getCommandString();
			options = parser.getOptionString();
										
			if (commandName != "") {
                numCommandsRun++;
                mout->setExecuting(true);
                mout->resetCommandErrors();
                
                //executes valid command
                mout->setChangedSeqNames(true);
               
                Command* command = cFactory->getCommand(commandName, options);
                quitCommandCalled = command->execute();
                delete command;
                
                //if we aborted command
                if (quitCommandCalled == 2) {  mout->mothurOut("[ERROR]: did not complete " + commandName + ".\n");  }
                
                if (mout->getControl_pressed()) { break;  }
                mout->setControl_pressed(false);
                mout->setExecuting(false);
                                
            }else {	mout->mothurOut("[ERROR]: Invalid.\n"); }
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
		bool ignoreSemiColons = false;
		
		//go through string until you reach ; or end
		while (count < commandString.length()) { 
			
			 //you want to ignore any ; until you reach the next '
			if ((commandString[count] == '\'') && (!ignoreSemiColons)) {  ignoreSemiColons = true;  } 
			else if ((commandString[count] == '\'') && (ignoreSemiColons)) {  ignoreSemiColons = false;  } 
				
			if ((commandString[count] == ';') && (!ignoreSemiColons)) {  break;   }
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
