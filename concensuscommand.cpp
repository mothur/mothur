/*
 *  concensuscommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 4/29/09.
 *  Copyright 2009 Schloss Lab UMASS AMherst. All rights reserved.
 *
 */

#include "concensuscommand.h"

//**********************************************************************************************************************

ConcensusCommand::ConcensusCommand(){
	try {
		globaldata = GlobalData::getInstance();
		t = globaldata->gTree;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ConcensusCommand class Function ConcensusCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ConcensusCommand class function ConcensusCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}
//**********************************************************************************************************************

ConcensusCommand::~ConcensusCommand(){}

//**********************************************************************************************************************

int ConcensusCommand::execute(){
	try {
		
		if (t.size() == 0) { return 0; }
		else {  
			numNodes = t[0]->getNumNodes();
			numLeaves = t[0]->getNumLeaves();
		}
		
		//initialize nodepairs
		nodePairs.resize(numNodes-numLeaves);
		
		//process each tree and fill counts
		for (int i = 0; i < t.size(); i++) {
			//process each nonleaf node
			for (int j = numLeaves; j < numNodes; j++) {
				//get pairing
				int leftChild = t[i]->tree[j].getLChild();
				int rightChild = t[i]->tree[j].getRChild();
				string pair = toString(leftChild) + "-" + toString(rightChild);
				
				//if this is an existing pairing for this node then increment the count otherwise add new pairing.
				it = nodePairs[j-numLeaves].find(pair);
				if (it != nodePairs[j-numLeaves].end()) {//already have that score
					nodePairs[j-numLeaves][pair]++;
				}else{//first time we have seen this score
					nodePairs[j-numLeaves][pair] = 1;
				}
			}
		}
		
		
		//print out pairings for testing
		/*for (int i = 0; i < nodePairs.size(); i++) {
			cout << "pairs for node " << i+numLeaves << endl;
			for (it = nodePairs[i].begin(); it != nodePairs[i].end(); it++) {
				cout << " pair = " << it->first <<  " count = " << it->second << endl;
			}
		}*/
		
		//open file for pairing not included in the tree
		notIncluded = getRootName(globaldata->inputFileName) + "concensuspairs";
		openOutputFile(notIncluded, out2);
		
		concensusTree = new Tree();
		
		//set relationships for nonleaf nodes
		for (int j = numLeaves; j < numNodes; j++) {
			
			//find that nodes pairing with the highest count
			int large = 0;
			for (it = nodePairs[j-numLeaves].begin(); it != nodePairs[j-numLeaves].end(); it++) {
				if (it->second > large) { large = it->second;  it2 = it; }
			}
			
			string pair = it2->first;
			int pos = pair.find_first_of('-');
			string lefts =  pair.substr(0, pos);
			string rights =  pair.substr(pos+1, pair.length());

			//converts string to int
			int left = atoi(lefts.c_str());
			int right = atoi(rights.c_str());
			
			//set parents and children
			concensusTree->tree[j].setChildren(left, right);
			concensusTree->tree[left].setParent(j);
			concensusTree->tree[right].setParent(j);
			
			// set Branchlengths //
			//if your children are leaves remove their branchlengths
			if (concensusTree->tree[left].getLChild() == -1) {  concensusTree->tree[left].setBranchLength(1.0); }
			if (concensusTree->tree[right].getLChild() == -1) {  concensusTree->tree[right].setBranchLength(1.0); }
			
			//set your branch length to the percentage of trees with this pairing
			float bLength = (float) it2->second / (float) t.size();
			concensusTree->tree[j].setBranchLength(bLength);
			
			//print out set used
			string leftName, rightName;
			getNames(it2->first, leftName, rightName);
			
			out2 << "Node " << j+1 << " in concensus tree: " << leftName << "-" << rightName << '\t' << (float)it2->second / (float) t.size() << endl; 
			out2 << "Node " << j+1 << " sets not included in concensus tree: " << endl; 
			
			//print out sets not included
			for (it = nodePairs[j-numLeaves].begin(); it != nodePairs[j-numLeaves].end(); it++) {
				if (it != it2) { 
					getNames(it->first, leftName, rightName);
					out2 << leftName << "-" << rightName << '\t' << (float)it->second / (float) t.size() << endl; 
				}
			}
			out2 << endl;

		}
		
		concensusTree->assembleTree();
		
		outputFile = getRootName(globaldata->inputFileName) + "concensus.tre";
		openOutputFile(outputFile, out);
		
		concensusTree->print(out);
		
		out.close(); out2.close();
		
		delete concensusTree;
		
		return 0;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ConcensusCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ConcensusCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		
}
//**********************************************************************************************************************

void ConcensusCommand::getNames(string pair, string& leftName, string& rightName) {
	try {
		int pos = pair.find_first_of('-');
		string lefts =  pair.substr(0, pos);
		string rights =  pair.substr(pos+1, pair.length());

		//converts string to int
		int	left = atoi(lefts.c_str());
		int	right = atoi(rights.c_str());
					
		//get name
		leftName = concensusTree->tree[left].getName();  
		rightName = concensusTree->tree[right].getName();
		 
		//if you are not a leaf use node number as name
		if (leftName == "") {  leftName = toString(left+1); }
		if (rightName == "") {  rightName = toString(right+1); }
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ConcensusCommand class Function getNames. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ConcensusCommand class function getNames. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		
}

//**********************************************************************************************************************


