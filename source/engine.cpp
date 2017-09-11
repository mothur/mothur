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

/***********************************************************************/

InteractEngine::InteractEngine(string path){

	
	string temppath = path.substr(0, (path.find_last_of("othur")-5));
	
	//this will happen if you set the path variable to contain mothur's exe location
	if (temppath == "") { path = mout->findProgramPath("mothur"); }
    else { path = temppath; }
	
	mout->setProgramPath(mout->getFullPathName(path));
    mout->setBlastPath(mout->getProgramPath());

    //if you haven't set your own location
    #ifdef MOTHUR_FILES
    #else
        //set default location to search for files to mothur's executable location.  This will resolve issue of double-clicking on the executable which opens mothur and sets pwd to your home directory instead of the mothur directory and leads to "unable to find file" errors.
        if (mout->getProgramPath() != "") { mout->setDefaultPath(mout->getProgramPath()); }
    #endif
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

                    
			if (mout->getChangedSeqNames()) { mout->mothurOut("[WARNING]: your sequence names contained ':'.  I changed them to '_' to avoid problems in your downstream analysis.\n"); }
                    
			input = getCommand();	
			
			if (mout->getControl_pressed()) { input = "quit()"; }
			
			//allow user to omit the () on the quit command
			if (input == "quit") { input = "quit()"; }

			CommandOptionParser parser(input);
			commandName = parser.getCommandString();
	
			options = parser.getOptionString();
			
			if (commandName != "") {
					mout->setExecuting(true);

					//executes valid command
                    mout->setChangedSeqNames(false);
					mout->setRunParse(true);
					mout->clearGroups();
					mout->clearAllGroups();
                    vector<string> temp;
					mout->setTreenames(temp);
					mout->setSaveNextLabel("");
                    mout->setCommandInputsConvertError(false);
					mout->setPrintedSharedHeaders(false);
					mout->setCurrentSharedBinLabels(temp);
					mout->setSharedBinLabelsInFile(temp);
                    mout->setPrintedListHeaders(false);
                    mout->setListBinLabelsInFile(temp);
							
					Command* command = cFactory->getCommand(commandName, options);
					if (mout->getCommandInputsConvertError()) { quitCommandCalled = 2; }
					else { quitCommandCalled = command->execute(); }
							
					//if we aborted command
					if (quitCommandCalled == 2) {  mout->mothurOut("[ERROR]: did not complete " + commandName + ".\n");  }

					mout->setControl_pressed(false);
					mout->setExecuting(false);
										
				}else {		
					mout->mothurOut("Invalid.\n");
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
	
		#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
			#ifdef USE_READLINE
				char* nextCommand = NULL;
				nextCommand = readline("\nmothur > ");
				
				if(nextCommand != NULL) {  add_history(nextCommand);  }	
				else{ //^D causes null string and we want it to quit mothur
					nextCommand = strdup("quit");
					mout->mothurOut(nextCommand);
                    mout->mothurOut("\n");
				}	
				
				mout->mothurOutJustToLog("\nmothur > " + toString(nextCommand) + "\n");
				return nextCommand;
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
		mout->errorOut(e, "Engine", "getCommand");
		exit(1);
	}
}
/***********************************************************************/
//This function opens the batchfile to be used by BatchEngine::getInput.
BatchEngine::BatchEngine(string path, string batchFileName){
	try {
	
		openedBatch = mout->openInputFile(batchFileName, inputBatchFile);
		
		string temppath = path.substr(0, (path.find_last_of("othur")-5));
	
		//this will happen if you set the path variable to contain mothur's exe location
		if (temppath == "") { path = mout->findProgramPath("mothur"); }
        else { path = temppath; }
		
        mout->setProgramPath(mout->getFullPathName(path));
        mout->setBlastPath(mout->getProgramPath());
        
        //if you haven't set your own location
#ifdef MOTHUR_FILES
#else
        //set default location to search for files to mothur's executable location.  This will resolve issue of double-clicking on the executable which opens mothur and sets pwd to your home directory instead of the mothur directory and leads to "unable to find file" errors.
        if (mout->getProgramPath() != "") { mout->setDefaultPath(mout->getProgramPath()); }
#endif

				
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
			mout->mothurOut("unable to open batchfile\n");
			return 1; 
		}
	
		string input = "";
		string commandName = "";
		string options = "";
		
		//CommandFactory cFactory;
		int quitCommandCalled = 0;
	    int count = 0;
		while(quitCommandCalled != 1){
			
			input = getNextCommand(inputBatchFile);
			count++;
			
            if (input[0] != '#') {
				if (mout->getChangedSeqNames()) { mout->mothurOut("[WARNING]: your sequence names contained ':'.  I changed them to '_' to avoid problems in your downstream analysis.\n"); }
				
				mout->mothurOut("\nmothur > " + input + "\n");
				
							
				if (mout->getControl_pressed()) { input = "quit()"; }
				
				//allow user to omit the () on the quit command
				if (input == "quit") { input = "quit()"; }

				CommandOptionParser parser(input);
				commandName = parser.getCommandString();
				options = parser.getOptionString();
										
				if (commandName != "") {
					mout->setExecuting(true);
					
					//executes valid command
                    mout->setChangedSeqNames(false);
					mout->setRunParse(true);
					mout->clearGroups();
					mout->clearAllGroups();
					vector<string> temp;
					mout->setTreenames(temp);
					mout->setSaveNextLabel("");
					mout->setCommandInputsConvertError(false);
                    mout->setPrintedSharedHeaders(false);
                    mout->setPrintedListHeaders(false);
                    mout->setListBinLabelsInFile(temp);
                    mout->setCurrentSharedBinLabels(temp);
                    mout->setSharedBinLabelsInFile(temp);
							
					Command* command = cFactory->getCommand(commandName, options);
					if (mout->getCommandInputsConvertError()) { quitCommandCalled = 2; }
					else { quitCommandCalled = command->execute(); }
							
					//if we aborted command
					if (quitCommandCalled == 2) {  mout->mothurOut("[ERROR]: did not complete " + commandName + ".\n");  }
                    
                    if (mout->getControl_pressed()) { break;  }
                    mout->setControl_pressed(false);
                    mout->setExecuting(false);
										
				}else {		
					mout->mothurOut("Invalid.\n");
				}
				
			}
			mout->gobble(inputBatchFile);
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
		else { nextcommand = mout->getline(inputBatchFile); }
		
		return nextcommand;
	}
	catch(exception& e) {
		mout->errorOut(e, "BatchEngine", "getNextCommand");
		exit(1);
	}
}

/***********************************************************************/
/***********************************************************************/
//This function opens the batchfile to be used by BatchEngine::getInput.
ScriptEngine::ScriptEngine(string path, string commandString){
	try {
		
		//remove quotes
		listOfCommands = commandString.substr(1, (commandString.length()-1));
		
		string temppath = path.substr(0, (path.find_last_of("othur")-5));

		//this will happen if you set the path variable to contain mothur's exe location
		if (temppath == "") { path = mout->findProgramPath("mothur"); }
        else { path = temppath; }
		
        mout->setProgramPath(mout->getFullPathName(path));
        mout->setBlastPath(mout->getProgramPath());
    
        //if you haven't set your own location
#ifdef MOTHUR_FILES
#else
        //set default location to search for files to mothur's executable location.  This will resolve issue of double-clicking on the executable which opens mothur and sets pwd to your home directory instead of the mothur directory and leads to "unable to find file" errors.
        //string tempDefault = path.substr(0, (path.find_last_of('m')));
        if (path != "") { mout->setDefaultPath(path); }
#endif

				
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
	
		while(quitCommandCalled != 1){
			
			input = getNextCommand(listOfCommands);	
			
			if (input == "") { input = "quit()"; }
                    
            if (mout->getChangedSeqNames()) { mout->mothurOut("[WARNING]: your sequence names contained ':'.  I changed them to '_' to avoid problems in your downstream analysis.\n"); }
			
			if (mout->getGui()) {
				if ((input.find("quit") != string::npos) || (input.find("set.logfile") != string::npos)) {}
				else if ((input.find("get.current") != string::npos) && (!mout->hasCurrentFiles())) {}
				else {  mout->mothurOut("\nmothur > " + input + "\n");  }
			}else{
				mout->mothurOut("\nmothur > " + input + "\n");
			}
			
            if (mout->getControl_pressed()) { input = "quit()"; }
				
			//allow user to omit the () on the quit command
			if (input == "quit") { input = "quit()"; }

			CommandOptionParser parser(input);
			commandName = parser.getCommandString();
			options = parser.getOptionString();
										
			if (commandName != "") {
					mout->setExecuting(true);
					
					//executes valid command
                    mout->setChangedSeqNames(false);
					mout->setRunParse(true);
					mout->clearGroups();
					mout->clearAllGroups();
					vector<string> temp;
					mout->setTreenames(temp);
					mout->setSaveNextLabel("");
                    mout->setCommandInputsConvertError(false);
                    mout->setPrintedSharedHeaders(false);
                    mout->setCurrentSharedBinLabels(temp);
                    mout->setSharedBinLabelsInFile(temp);
                    mout->setPrintedListHeaders(false);
                    mout->setListBinLabelsInFile(temp);

					Command* command = cFactory->getCommand(commandName, options);
					if (mout->getCommandInputsConvertError()) { quitCommandCalled = 2; }
					else { quitCommandCalled = command->execute(); }
					
					//if we aborted command
					if (quitCommandCalled == 2) {  mout->mothurOut("[ERROR]: did not complete " + commandName + ".\n");  }
					
                    if (mout->getControl_pressed()) { break;  }
                    mout->setControl_pressed(false);
                    mout->setExecuting(false);
									
				}else {		
					mout->mothurOut("Invalid.\n");
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
