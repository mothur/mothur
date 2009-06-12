/*
 *  GetlineCommand.cpp
 *  Mothur
 *
 *  Created by Thomas Ryabin on 1/30/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "getlinecommand.h"

//**********************************************************************************************************************
GetlineCommand::GetlineCommand(string option){
	try {
		globaldata = GlobalData::getInstance();
		abort = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			if (option != "") { cout << "There are no valid parameters for the get.line command." << endl; abort = true; }
			
			if ((globaldata->getListFile() == "") && (globaldata->getRabundFile() == "") && (globaldata->getSabundFile() == "")) { cout << "You must read a list, sabund or rabund before you can use the get.line command." << endl; abort = true; }				
		}

	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the GetlineCommand class Function GetlineCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the GetlineCommand class function GetlineCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
			
}
//**********************************************************************************************************************

void GetlineCommand::help(){
	try {
		cout << "The get.line command can only be executed after a successful read.otu command." << "\n";
		cout << "You may not use any parameters with the get.line command." << "\n";
		cout << "The get.line command should be in the following format: " << "\n";
		cout << "get.line()" << "\n";
		cout << "Example get.line()." << "\n";
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the GetlineCommand class Function help. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the GetlineCommand class function help. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

//**********************************************************************************************************************

GetlineCommand::~GetlineCommand(){
}

//**********************************************************************************************************************

int GetlineCommand::execute(){
	try {
	
		if (abort == true) { return 0; }
	
		filename = globaldata->inputFileName;
		ifstream in;
		openInputFile(filename, in);
		string label;
		int numBins = 0;
		int count = -1;
		int line = 1;
		while(in.good()) {
			if(count > numBins)
				count = 0;
			if(count == 0) {
				cout << line << "\n";
				in >> numBins;
				line++;
			}
			in >> label;
			count++;
		}
		return 0;		
	}

	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the GetlineCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the GetlineCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}



