/*
 *  readtreecommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 1/23/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "readtreecommand.h"

//**********************************************************************************************************************
ReadTreeCommand::ReadTreeCommand(){
	try {
		globaldata = GlobalData::getInstance();
		filename = globaldata->inputFileName;
		
		//read in group map info.
		treeMap = new TreeMap(globaldata->getGroupFile());
		treeMap->readMap();
		globaldata->gTreemap = treeMap;

		read = new ReadNewickTree(filename);
		
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ReadTreeCommand class Function ReadTreeCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ReadTreeCommand class function ReadTreeCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

//**********************************************************************************************************************

ReadTreeCommand::~ReadTreeCommand(){
	delete read;
}

//**********************************************************************************************************************

int ReadTreeCommand::execute(){
	try {
	
		read->read(); 
		
		vector<Tree*> T = globaldata->gTree;
		
		//assemble users trees
		for (int i = 0; i < T.size(); i++) {
			T[i]->assembleTree();
		}
		
		return 0;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ReadTreeCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ReadTreeCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

//**********************************************************************************************************************
