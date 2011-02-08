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
	
	string temppath = path.substr(0, (path.find_last_of('m')));
	
	//this will happen if you set the path variable to contain mothur's exe location
	if (temppath == "") { 
	
		string envPath = getenv("PATH");
		
		//delimiting path char
		char delim;
		#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
			delim = ':';
		#else
			delim = ';';
		#endif
		
		//break apart path variable by ':'
		vector<string> dirs;
		mout->splitAtChar(envPath, dirs, delim);
		
		//get path related to mothur
		string mothurPath = "";
		for (int i = 0; i < dirs.size(); i++) {
			//to lower so we can find it
			string tempLower = "";
			for (int j = 0; j < dirs[i].length(); j++) {  tempLower += tolower(dirs[i][j]);  }
			
			//is this mothurs path?
			if (tempLower.find("mothur") != -1) {  mothurPath = dirs[i]; break;  }
		}
		
		//add mothur so it looks like what argv would look like
		#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
			mothurPath += "/mothur";
		#else
			mothurPath += "\\mothur";
		#endif
		
		path = mothurPath;
	}
	
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

			#ifdef USE_MPI
				int pid, processors;
				MPI_Status status;
				MPI_Comm_rank(MPI_COMM_WORLD, &pid); 
				MPI_Comm_size(MPI_COMM_WORLD, &processors);
				
				if (pid == 0) {
				
			#endif
			
			mout->mothurOutEndLine();
			
			input = getCommand();	
			mout->mothurOutEndLine();	
			
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
					Command* command = cFactory->getCommand(commandName, options);
					quitCommandCalled = command->execute();
							
					//if we aborted command
					if (quitCommandCalled == 2) {  mout->mothurOut("[ERROR]: did not complete " + commandName + "."); mout->mothurOutEndLine(); }

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
					strcpy(nextCommand, "quit"); 
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
	
		openedBatch = mout->openInputFile(batchFileName, inputBatchFile);
		
		string temppath = path.substr(0, (path.find_last_of('m')));
	
		//this will happen if you set the path variable to contain mothur's exe location
		if (temppath == "") { 
		
			string envPath = getenv("PATH");
			
			//delimiting path char
			char delim;
			#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
				delim = ':';
			#else
				delim = ';';
			#endif
			
			//break apart path variable by ':'
			vector<string> dirs;
			mout->splitAtChar(envPath, dirs, delim);
			
			//get path related to mothur
			string mothurPath = "";
			for (int i = 0; i < dirs.size(); i++) {
				//to lower so we can find it
				string tempLower = "";
				for (int j = 0; j < dirs[i].length(); j++) {  tempLower += tolower(dirs[i][j]);  }
				
				//is this mothurs path?
				if (tempLower.find("mothur") != -1) {  mothurPath = dirs[i]; break;  }
			}
			
			//add mothur so it looks like what argv would look like
			#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
				mothurPath += "/mothur";
			#else
				mothurPath += "\\mothur";
			#endif
			
			path = mothurPath;
		}

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
						
//cout << pid << " is here " << commandName << '\t' << count << endl;
						if ((cFactory->MPIEnabled(commandName)) || (pid == 0)) {
					#endif
					//executes valid command
					Command* command = cFactory->getCommand(commandName, options);
					quitCommandCalled = command->execute();
							
					//if we aborted command
					if (quitCommandCalled == 2) {  mout->mothurOut("[ERROR]: did not complete " + commandName + "."); mout->mothurOutEndLine(); }

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
		globaldata = GlobalData::getInstance();
		
		//remove quotes
		listOfCommands = commandString.substr(1, (commandString.length()-1));
		
		string temppath = path.substr(0, (path.find_last_of('m')));
	
		//this will happen if you set the path variable to contain mothur's exe location
		if (temppath == "") { 
		
			string envPath = getenv("PATH");
			
			//delimiting path char
			char delim;
			#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
				delim = ':';
			#else
				delim = ';';
			#endif
			
			//break apart path variable by ':'
			vector<string> dirs;
			mout->splitAtChar(envPath, dirs, delim);
			
			//get path related to mothur
			string mothurPath = "";
			for (int i = 0; i < dirs.size(); i++) {
				//to lower so we can find it
				string tempLower = "";
				for (int j = 0; j < dirs[i].length(); j++) {  tempLower += tolower(dirs[i][j]);  }
				
				//is this mothurs path?
				if (tempLower.find("mothur") != -1) {  mothurPath = dirs[i]; break;  }
			}
			
			//add mothur so it looks like what argv would look like
			#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
				mothurPath += "/mothur";
			#else
				mothurPath += "\\mothur";
			#endif
			
			path = mothurPath;
		}

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
	
		while(quitCommandCalled != 1){
			
			#ifdef USE_MPI
				int pid, processors;
				MPI_Status status;
				MPI_Comm_rank(MPI_COMM_WORLD, &pid); 
				MPI_Comm_size(MPI_COMM_WORLD, &processors);
				
				if (pid == 0) {
				
			#endif
			
			input = getNextCommand(listOfCommands);	
			
			if (input == "") { input = "quit()"; }
			
			mout->mothurOutEndLine();
			mout->mothurOut("mothur > " + input);
			mout->mothurOutEndLine();
			
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
					Command* command = cFactory->getCommand(commandName, options);
					quitCommandCalled = command->execute();
					
					//if we aborted command
					if (quitCommandCalled == 2) {  mout->mothurOut("[ERROR]: did not complete " + commandName + "."); mout->mothurOutEndLine(); }
							
					mout->control_pressed = 0;
					mout->executing = false;
									
					#ifdef USE_MPI
					//cout << pid << " is done in execute" << endl;
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
