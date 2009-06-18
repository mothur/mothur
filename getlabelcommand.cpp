/*
 *  GetlabelCommand.cpp
 *  Mothur
 *
 *  Created by Thomas Ryabin on 1/30/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "getlabelcommand.h"

//**********************************************************************************************************************

GetlabelCommand::GetlabelCommand(string option){
	try {
		globaldata = GlobalData::getInstance();
		abort = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			if (option != "") { cout << "There are no valid parameters for the get.label command." << endl; abort = true; }
			
			if ((globaldata->getListFile() == "") && (globaldata->getRabundFile() == "") && (globaldata->getSabundFile() == "")) { cout << "You must read a list, sabund or rabund before you can use the get.label command." << endl; abort = true; }				
		}

	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the GetlabelCommand class Function GetlabelCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the GetlabelCommand class function GetlabelCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
			
}
//**********************************************************************************************************************

void GetlabelCommand::help(){
	try {
		cout << "The get.label command can only be executed after a successful read.otu command." << "\n";
		cout << "You may not use any parameters with the get.label command." << "\n";
		cout << "The get.label command should be in the following format: " << "\n";
		cout << "get.label()" << "\n";
		cout << "Example get.label()." << "\n";
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the GetlabelCommand class Function help. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the GetlabelCommand class function help. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

//**********************************************************************************************************************

GetlabelCommand::~GetlabelCommand(){
}

//**********************************************************************************************************************

int GetlabelCommand::execute(){
	try {
		
		if (abort == true) { return 0; }
		
		filename = globaldata->inputFileName;
		ifstream in;
		openInputFile(filename, in);
		string label;
		int numBins = 0;
		int count = -1;
		while(in.good()) {
			if(count > numBins)
				count = 0;
			if(count == 0) {
				cout << label << "\n";
				in >> numBins;
			}
			in >> label;
			count++;
		}	
		
		in.close();
		return 0;	
	}

	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the GetlabelCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the GetlabelCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

