/*
 *  readsabundfilecommand.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "readsabundfilecommand.h"

//**********************************************************************************************************************
ReadSAbundFileCommand::ReadSAbundFileCommand(){
	try {
		globaldata = GlobalData::getInstance();
		filename = globaldata->inputFileName;	
		read = new ReadPhilFile(filename);
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ReadSAbundFileCommand class Function ReadSAbundFileCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ReadSAbundFileCommand class function ReadSAbundFileCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

//**********************************************************************************************************************

ReadSAbundFileCommand::~ReadSAbundFileCommand(){
	delete read;
}

//**********************************************************************************************************************

int ReadSAbundFileCommand::execute(){
	try {
		read->read(&*globaldata); 
		return 0;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ReadSAbundFileCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ReadSAbundFileCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
//**********************************************************************************************************************
