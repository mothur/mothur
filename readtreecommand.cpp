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
ReadTreeCommand::ReadTreeCommand(string option){
	try {
		globaldata = GlobalData::getInstance();
		abort = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"tree","group"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string, string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
		
			//check to make sure all parameters are valid for command
			for (map<string, string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			globaldata->newRead();
			
			//check for required parameters
			treefile = validParameter.validFile(parameters, "tree", true);
			if (treefile == "not open") { abort = true; }
			else if (treefile == "not found") { treefile = ""; cout << "tree is a required parameter for the read.tree command." << endl; abort = true;  }	
			else {  globaldata->setTreeFile(treefile);  globaldata->setFormat("tree"); 	}
			
			groupfile = validParameter.validFile(parameters, "group", true);
			if (groupfile == "not open") { abort = true; }	
			else if (groupfile == "not found") { groupfile = ""; cout << "group is a required parameter for the read.tree command." << endl; abort = true;	}
			else {  
				globaldata->setGroupFile(groupfile); 
				//read in group map info.
				treeMap = new TreeMap(groupfile);
				treeMap->readMap();
				globaldata->gTreemap = treeMap;
			}
			
			if (abort == false) {
				filename = treefile;
				read = new ReadNewickTree(filename);
			}
						
		}
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

void ReadTreeCommand::help(){
	try {
		cout << "The read.tree command must be run before you execute a unifrac.weighted, unifrac.unweighted. " << "\n";
		cout << "It also must be run before using the parsimony command, unless you are using the randomtree parameter." << "\n";
		cout << "The read.tree command should be in the following format: read.tree(tree=yourTreeFile, group=yourGroupFile)." << "\n";
		cout << "The tree and group parameters are both required." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. tree), '=' and parameters (i.e.yourTreefile)." << "\n" << "\n";
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ReadTreeCommand class Function help. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ReadTreeCommand class function help. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
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
	
		if (abort == true) { return 0; }
		
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
