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
	
	mout->argv = path;

    //if you haven't set your own location
    #ifdef MOTHUR_FILES
    #else
        //set default location to search for files to mothur's executable location.  This will resolve issue of double-clicking on the executable which opens mothur and sets pwd to your home directory instead of the mothur directory and leads to "unable to find file" errors.
        string tempDefault = path.substr(0, (path.find_last_of('m')));
        if (tempDefault != "") { mout->setDefaultPath(tempDefault); }
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

			#ifdef USE_MPI
				int pid, processors;
				MPI_Status status;
				MPI_Comm_rank(MPI_COMM_WORLD, &pid); 
				MPI_Comm_size(MPI_COMM_WORLD, &processors);
				
				if (pid == 0) {
				
			#endif
                    
			if (mout->changedSeqNames) { mout->mothurOut("[WARNING]: your sequence names contained ':'.  I changed them to '_' to avoid problems in your downstream analysis.\n"); }
                    
			input = getCommand();	
			
			if (mout->control_pressed) { input = "quit()"; }
			
			//allow user to omit the () on the quit command
			if (input == "quit") { input = "quit()"; }

			
			#ifdef USE_MPI
				//send commandName
				for(int i = 1; i < processors; i++) { 
						int length = input.length();
						MPI_Send(&length, 1, MPI_INT, i, 2001, MPI_COMM_WORLD);
						MPI_Send(&input[0], length, MPI_CHAR, i, 2001, MPI_COMM_WORLD);
	
					}
				}else {
					int length;
					MPI_Recv(&length, 1, MPI_INT, 0, 2001, MPI_COMM_WORLD, &status);
					//recieve container
					char* tempBuf = new char[length];
					MPI_Recv(&tempBuf[0], length, MPI_CHAR, 0, 2001, MPI_COMM_WORLD, &status);
					
					input = tempBuf;
					if (input.length() > length) { input = input.substr(0, length);  }
					delete tempBuf;	
				}

			
			#endif
		
			CommandOptionParser parser(input);
			commandName = parser.getCommandString();
	
			options = parser.getOptionString();
			
			if (commandName != "") {
					mout->executing = true;
					#ifdef USE_MPI
						int pid;
						MPI_Comm_rank(MPI_COMM_WORLD, &pid); 
						
						if ((cFactory->MPIEnabled(commandName)) || (pid == 0)) {
					//cout << pid << " is in execute " << commandName << endl;
					#endif
					//executes valid command
                    mout->changedSeqNames = false;
					mout->runParse = true;
					mout->clearGroups();
					mout->clearAllGroups();
					mout->Treenames.clear();
					mout->saveNextLabel = "";
                    mout->commandInputsConvertError = false;
					mout->printedSharedHeaders = false;
					mout->currentSharedBinLabels.clear();
					mout->sharedBinLabelsInFile.clear();
                    mout->printedListHeaders = false;
                    mout->listBinLabelsInFile.clear();
							
					Command* command = cFactory->getCommand(commandName, options);
					if (mout->commandInputsConvertError) { quitCommandCalled = 2; }
					else { quitCommandCalled = command->execute(); }
							
					//if we aborted command
					if (quitCommandCalled == 2) {  mout->mothurOut("[ERROR]: did not complete " + commandName + ".\n");  }

					mout->control_pressed = 0;
					mout->executing = false;
										
					#ifdef USE_MPI
						}
                        MPI_Barrier(MPI_COMM_WORLD); //make everyone wait - just in case
					#endif
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
		
		mout->argv = path;
        
        //if you haven't set your own location
#ifdef MOTHUR_FILES
#else
        //set default location to search for files to mothur's executable location.  This will resolve issue of double-clicking on the executable which opens mothur and sets pwd to your home directory instead of the mothur directory and leads to "unable to find file" errors.
        string tempDefault = path.substr(0, (path.find_last_of('m')));
        if (tempDefault != "") { mout->setDefaultPath(tempDefault); }
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
			
			#ifdef USE_MPI
				int pid, processors;
				MPI_Status status;
				MPI_Comm_rank(MPI_COMM_WORLD, &pid); 
				MPI_Comm_size(MPI_COMM_WORLD, &processors);
				
				if (pid == 0) {
				
			#endif
			
			input = getNextCommand(inputBatchFile);
			count++;
			
			#ifdef USE_MPI
				//send commandName
				for(int i = 1; i < processors; i++) { 
						int length = input.length();
						MPI_Send(&length, 1, MPI_INT, i, 2001, MPI_COMM_WORLD);
						MPI_Send(&input[0], length, MPI_CHAR, i, 2001, MPI_COMM_WORLD);
	
					}
				}else {
					int length;
					MPI_Recv(&length, 1, MPI_INT, 0, 2001, MPI_COMM_WORLD, &status);
					//recieve container
					char* tempBuf = new char[length];
					MPI_Recv(&tempBuf[0], length, MPI_CHAR, 0, 2001, MPI_COMM_WORLD, &status);
					
					input = tempBuf;
					if (input.length() > length) { input = input.substr(0, length);  }
					delete tempBuf;	
				}

			
			#endif

			
			
			if (input[0] != '#') {
				if (mout->changedSeqNames) { mout->mothurOut("[WARNING]: your sequence names contained ':'.  I changed them to '_' to avoid problems in your downstream analysis.\n"); }
				
				mout->mothurOut("\nmothur > " + input + "\n");
				
							
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
						
//cout << pid << " is here " << commandName << '\t' << count << endl;
						if ((cFactory->MPIEnabled(commandName)) || (pid == 0)) {
					#endif
					//executes valid command
                    mout->changedSeqNames = false;
					mout->runParse = true;
					mout->clearGroups();
					mout->clearAllGroups();
					mout->Treenames.clear();
					mout->saveNextLabel = "";
					mout->commandInputsConvertError = false;
                    mout->printedSharedHeaders = false;
                    mout->currentSharedBinLabels.clear();
                    mout->sharedBinLabelsInFile.clear();
                    mout->printedListHeaders = false;
                    mout->listBinLabelsInFile.clear();

							
					Command* command = cFactory->getCommand(commandName, options);
					if (mout->commandInputsConvertError) { quitCommandCalled = 2; }
					else { quitCommandCalled = command->execute(); }
							
					//if we aborted command
					if (quitCommandCalled == 2) {  mout->mothurOut("[ERROR]: did not complete " + commandName + ".\n");  }

					mout->control_pressed = 0;
					mout->executing = false;
										
					#ifdef USE_MPI
						}
                        MPI_Barrier(MPI_COMM_WORLD); //make everyone wait - just in case
					#endif
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
		
		mout->argv = path;
    
        //if you haven't set your own location
#ifdef MOTHUR_FILES
#else
        //set default location to search for files to mothur's executable location.  This will resolve issue of double-clicking on the executable which opens mothur and sets pwd to your home directory instead of the mothur directory and leads to "unable to find file" errors.
        string tempDefault = path.substr(0, (path.find_last_of('m')));
        if (tempDefault != "") { mout->setDefaultPath(tempDefault); }
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
			
			#ifdef USE_MPI
				int pid, processors;
				MPI_Status status;
				MPI_Comm_rank(MPI_COMM_WORLD, &pid); 
				MPI_Comm_size(MPI_COMM_WORLD, &processors);
				
				if (pid == 0) {
				//cout << pid << " is here " << processors << endl;
			#endif
			
			input = getNextCommand(listOfCommands);	
			
			if (input == "") { input = "quit()"; }
                    
            if (mout->changedSeqNames) { mout->mothurOut("[WARNING]: your sequence names contained ':'.  I changed them to '_' to avoid problems in your downstream analysis.\n"); }
			
			if (mout->gui) {
				if ((input.find("quit") != string::npos) || (input.find("set.logfile") != string::npos)) {}
				else if ((input.find("get.current") != string::npos) && (!mout->hasCurrentFiles())) {}
				else {  mout->mothurOut("\nmothur > " + input + "\n");  }
			}else{
				mout->mothurOut("\nmothur > " + input + "\n");
			}
			
			#ifdef USE_MPI
				//send commandName
				for(int i = 1; i < processors; i++) {
                    //cout << pid << " is here " << input << endl;
						int length = input.length();
						MPI_Send(&length, 1, MPI_INT, i, 2001, MPI_COMM_WORLD);
                    //cout << pid << " is here " << length << '\t' << input << endl;
						MPI_Send(&input[0], length, MPI_CHAR, i, 2001, MPI_COMM_WORLD);
	//cout << pid << " is here " << length << '\t' << input << endl;
					}
				}else {
					int length;
					MPI_Recv(&length, 1, MPI_INT, 0, 2001, MPI_COMM_WORLD, &status);
                    //cout << pid << " is here " << length << endl;
					//recieve container
					char* tempBuf = new char[length];
					MPI_Recv(&tempBuf[0], length, MPI_CHAR, 0, 2001, MPI_COMM_WORLD, &status);
					//cout << pid << " is here " << length << '\t' << tempBuf << endl;
					input = tempBuf;
					if (input.length() > length) { input = input.substr(0, length);  }
					delete tempBuf;	
				}

			
			#endif
			
			
			if (mout->control_pressed) { input = "quit()"; }
				
			//allow user to omit the () on the quit command
			if (input == "quit") { input = "quit()"; }

			CommandOptionParser parser(input);
			commandName = parser.getCommandString();
			options = parser.getOptionString();
										
			if (commandName != "") {
					mout->executing = true;
					#ifdef USE_MPI
						int pid, numProcesses;
						
						MPI_Comm_rank(MPI_COMM_WORLD, &pid); 
						MPI_Comm_size(MPI_COMM_WORLD, &numProcesses); 
					
//cout << pid << " is here " << commandName  << endl;
						if ((cFactory->MPIEnabled(commandName)) || (pid == 0)) {
							//cout << pid << " is in execute" << endl;
					#endif
					//executes valid command
                    mout->changedSeqNames = false;
					mout->runParse = true;
					mout->clearGroups();
					mout->clearAllGroups();
					mout->Treenames.clear();
					mout->saveNextLabel = "";
                    mout->commandInputsConvertError = false;
                    mout->printedSharedHeaders = false;
                    mout->currentSharedBinLabels.clear();
                    mout->sharedBinLabelsInFile.clear();
                    mout->printedListHeaders = false;
                    mout->listBinLabelsInFile.clear();

					Command* command = cFactory->getCommand(commandName, options);
					if (mout->commandInputsConvertError) { quitCommandCalled = 2; }
					else { quitCommandCalled = command->execute(); }
					
					//if we aborted command
					if (quitCommandCalled == 2) {  mout->mothurOut("[ERROR]: did not complete " + commandName + ".\n");  }
							
					mout->control_pressed = 0;
					mout->executing = false;
									
					#ifdef USE_MPI
					//cout << pid << " is done in execute" << endl;
						}
                        MPI_Barrier(MPI_COMM_WORLD); //make everyone wait - just in case
					#endif
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
