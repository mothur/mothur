/*
 *  readotu.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 1/20/09.
 *  Copyright 2009 Schloss Lab UMASS AMherst. All rights reserved.
 *
 */

#include "readotucommand.h"

//**********************************************************************************************************************
ReadOtuCommand::ReadOtuCommand(){
	try {
		globaldata = GlobalData::getInstance();
		filename = globaldata->inputFileName;
		read = new ReadOTUFile(filename);
		if (globaldata->getFormat() == "shared") {
			//read in group map info.
			groupMap = new GroupMap(globaldata->getGroupFile());
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ReadOtuCommand class Function ReadOtuCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ReadOtuCommand class function ReadOtuCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

//**********************************************************************************************************************

ReadOtuCommand::~ReadOtuCommand(){
	delete read;
}

//**********************************************************************************************************************

int ReadOtuCommand::execute(){
	try {
		read->read(&*globaldata); 
		if (globaldata->getFormat() == "shared") {
			groupMap->readMap();
			
			//if (globaldata->gGroupmap != NULL) { delete globaldata->gGroupmap;  }
			globaldata->gGroupmap = groupMap;		
			shared = new SharedCommand();
			shared->execute();

			parselist = new ParseListCommand();
			parselist->execute();
			
			//change format to shared  to speed up commands
			globaldata->setFormat("sharedfile");
			globaldata->setListFile("");
			globaldata->setGroupFile("");
			globaldata->setSharedFile(getRootName(filename) + "shared");
		}
		return 0;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ReadOtuCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ReadOtuCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
//**********************************************************************************************************************
