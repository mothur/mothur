/*
 *  sharedcommand.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "sharedcommand.h"

//**********************************************************************************************************************

SharedCommand::SharedCommand(){
	try {
		globaldata = GlobalData::getInstance();
		
		//getting output filename
		filename = globaldata->inputFileName;
		filename = getRootName(filename);
		filename = filename + "shared";
		openOutputFile(filename, out);
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SharedCommand class Function SharedCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SharedCommand class function SharedCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}

}
//**********************************************************************************************************************

int SharedCommand::execute(){
	try {
		globaldata = GlobalData::getInstance();
			
		//read in listfile
		read = new ReadPhilFile(globaldata->inputFileName);	
		read->read(&*globaldata); 
		input = globaldata->ginput;
		SharedList = globaldata->gSharedList;
		
		shared = new Shared();
		int i = 0;
		while(SharedList != NULL){
			shared->getSharedVectors(i, SharedList); //fills sharedGroups with new info and updates sharedVector
			SharedList = input->getSharedListVector(); //get new list vector to process
			printSharedData(); //prints info to the .shared file
			i++;
		}
		return 0;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SharedCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SharedCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}

}
//**********************************************************************************************************************
void SharedCommand::printSharedData() {
	try {
		//prints out horizontally
		for (it = shared->sharedGroups.begin(); it != shared->sharedGroups.end(); it++) {
			out << it->second->getLabel() << "\t" << it->first << "\t"; //prints out label and groupname
			it->second->print(out); // prints sharedrabundvector
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SharedCommand class Function printSharedData. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SharedCommand class function printSharedData. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	
}

//**********************************************************************************************************************

SharedCommand::~SharedCommand(){
	//delete list;
	delete read;
}

//**********************************************************************************************************************
