/*
 *  readsharedcommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 2/12/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "readsharedcommand.h"

//**********************************************************************************************************************
ReadSharedCommand::ReadSharedCommand(){
	try {
		globaldata = GlobalData::getInstance();
		filename = globaldata->inputFileName;
		read = new ReadPhilFile(filename);
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ReadSharedCommand class Function ReadSharedCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ReadSharedCommand class function ReadSharedCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

//**********************************************************************************************************************

ReadSharedCommand::~ReadSharedCommand(){
	delete read;
}

//**********************************************************************************************************************

int ReadSharedCommand::execute(){
	try {
		read->read(&*globaldata); 
	
		return 0;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ReadSharedCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ReadSharedCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
//**********************************************************************************************************************