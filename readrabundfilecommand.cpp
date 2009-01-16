/*
 *  readrabundfilecommand.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "readrabundfilecommand.h"

//**********************************************************************************************************************
ReadRAbundFileCommand::ReadRAbundFileCommand(){
	try {
		globaldata = GlobalData::getInstance();
		filename = globaldata->inputFileName;
		read = new ReadPhilFile(filename);
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ReadRAbundFileCommand class Function ReadRAbundFileCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ReadRAbundFileCommand class function ReadRAbundFileCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

//**********************************************************************************************************************

ReadRAbundFileCommand::~ReadRAbundFileCommand(){
	delete read;
}

//**********************************************************************************************************************

int ReadRAbundFileCommand::execute(){
	try {
		read->read(&*globaldata); 
		return 0;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ReadRAbundFileCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ReadRAbundFileCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
//**********************************************************************************************************************
