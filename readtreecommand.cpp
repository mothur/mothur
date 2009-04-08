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
		
		//memory leak prevention
		//if (globaldata->gTreemap != NULL) { delete globaldata->gTreemap;  }
		globaldata->gTreemap = treeMap;
		
		//get names in tree
		globaldata->parseTreeFile();

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
		int readOk;
		
		readOk = read->read(); 
		
		if (readOk != 0) { cout << "Read Terminated." << endl; globaldata->gTree.clear(); delete globaldata->gTreemap; return 0; }
		
		vector<Tree*> T = globaldata->gTree;
		
		//assemble users trees
		for (int i = 0; i < T.size(); i++) {
			T[i]->assembleTree();
		}

		//output any names that are in names file but not in tree
		if (globaldata->Treenames.size() < treeMap->getNumSeqs()) {
			for (int i = 0; i < treeMap->namesOfSeqs.size(); i++) {
				//is that name in the tree?
				int count = 0;
				for (int j = 0; j < globaldata->Treenames.size(); j++) {
					if (treeMap->namesOfSeqs[i] == globaldata->Treenames[j]) { break; } //found it
					count++;
				}
				
				//then you did not find it so report it 
				if (count == globaldata->Treenames.size()) { 
					cout << treeMap->namesOfSeqs[i] << " is in your namefile and not in your tree. It will be disregarded." << endl;
				}
			}
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
