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
	
	system("clear");
//	char buffer = ' ';
//	ifstream header("introtext.txt");
//	while(!header.eof()){
//		cout << buffer;
//		buffer = header.get();
//	}
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
		//bool errorFree;
		//ErrorCheck* errorCheckor = new ErrorCheck();
		
		cout << "mothur v.1.3.0" << endl;
		cout << "Last updated: 5/29/2009" << endl << endl;
		cout << "by" << endl;
		cout << "Patrick D. Schloss" << endl << endl;
		cout << "Department of Microbiology" << endl;
		cout << "The University of Massachusetts" << endl;
		cout << "pschloss@micro.umass.edu" << endl;
		cout << "http://schloss.micro.umass.edu/mothur" << endl << endl << endl;
		cout << "Distributed under the GNU General Public License" << endl << endl;
		cout << "Type 'help()' for information on the commands that are available" << endl << endl;
		cout << "Type 'quit()' to exit program" << endl;

		while(quitCommandCalled != 1){

			cout << endl << "mothur > ";
			getline(cin, input);
			if (cin.eof()) { input = "quit()"; }
			
			//allow user to omit the () on the quit command
			if (input == "quit") { input = "quit()"; }
			
			//errorFree = errorCheckor->checkInput(input);
			//if (errorFree == true) {
			CommandOptionParser parser(input);
			commandName = parser.getCommandString();
			options = parser.getOptionString();
			
			if (commandName != "") {
			
				//executes valid command
				CommandFactory cFactory;
				Command* command = cFactory.getCommand(commandName, options);
				quitCommandCalled = command->execute();
				
			}else {
				cout << "Your input contains errors. Please try again." << endl;
			}
		}	
		return 1;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the InteractEngine class Function getInput. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the InteractEngine class function getInput. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
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

		system("clear");
	
	//	char buffer = ' ';
	//	ifstream header("introtext.txt");
	//	while(!header.eof()){
	//		cout << buffer;
	//		buffer = header.get();
	//	}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the BatchEngine class Function BatchEngine. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the BatchEngine class function BatchEngine. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
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
		if (openedBatch == 1) {  cout << "unable to open batchfile" << endl;  return 1; }
	
		string input = "";
		string commandName = "";
		string options = "";
		
		//CommandFactory cFactory;
		int quitCommandCalled = 0;
	
		while(quitCommandCalled == 0){
		
			if (inputBatchFile.eof()) { input = "quit()"; }
			else { getline(inputBatchFile, input); }
			
			if (input[0] != '#') {
			
				cout << endl << "mothur > " << input << endl;
				
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
				}else {		cout << "Invalid." << endl;		}
				
			}
			gobble(inputBatchFile);
		}
		return 1;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the BatchEngine class Function getInput. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the BatchEngine class function getInput. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}


/***********************************************************************/

