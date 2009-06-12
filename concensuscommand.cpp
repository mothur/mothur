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

ConcensusCommand::ConcensusCommand(string option){
	try {
		globaldata = GlobalData::getInstance();
		abort = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			if (option != "") { cout << "There are no valid parameters for the concensus command." << endl; abort = true; }
			
			//no trees were read
			if (globaldata->gTree.size() == 0) {	cout << "You must execute the read.tree command, before you may use the concensus command." << endl; abort = true;  }
			else {  t = globaldata->gTree;	}
		}
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

void ConcensusCommand::help(){
	try {
		cout << "The concensus command can only be executed after a successful read.tree command." << "\n";
		cout << "The concensus command has no parameters." << "\n";
		cout << "The concensus command should be in the following format: concensus()." << "\n";
		cout << "The concensus command output two files: .concensus.tre and .concensuspairs." << "\n";
		cout << "The .concensus.tre file contains the concensus tree of the trees in your input file." << "\n";
		cout << "The branch lengths are the percentage of trees in your input file that had the given pair." << "\n";
		cout << "The .concensuspairs file contains a list of the internal nodes in your tree.  For each node, the pair that was used in the concensus tree " << "\n";
		cout << "is reported with its percentage, as well as the other pairs that were seen for that node but not used and their percentages." << "\n" << "\n";		
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ConcensusCommand class Function help. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ConcensusCommand class function help. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

//**********************************************************************************************************************

ConcensusCommand::~ConcensusCommand(){}

//**********************************************************************************************************************

int ConcensusCommand::execute(){
	try {
		
		if (abort == true) { return 0; }
		else {  
			numNodes = t[0]->getNumNodes();
			numLeaves = t[0]->getNumLeaves();
		}
		
		//get the possible pairings
		getSets();		
		
		//print out pairings for testing
		/*cout << "possible pairing " <<  endl;
		for (it2 = nodePairs.begin(); it2 != nodePairs.end(); it2++) {
			for (int i = 0; i < it2->first.size(); i++) {
				cout << it2->first[i] << " ";
			}
			cout << '\t' << it2->second << endl;
		}*/
		
		
		//open file for pairing not included in the tree
		notIncluded = getRootName(globaldata->inputFileName) + "concensuspairs";
		openOutputFile(notIncluded, out2);
		
		concensusTree = new Tree();
		
		it2 = nodePairs.find(treeSet);
		
		nodePairsInTree[treeSet] = it2->second; 
		
		//erase treeset because you are adding it
		nodePairs.erase(treeSet);
		
		//set count to numLeaves;
		count = numLeaves;
		
		buildConcensusTree(treeSet);
		
		concensusTree->assembleTree();
		
		//output species in order
		out2 << "Species in Order: " << endl << endl;
		for (int i = 0; i < treeSet.size(); i++) {  out2 << i+1 << ".  " << treeSet[i] << endl; }
		
		vector<string> temp; 
		//output sets included
		out2 << endl << "Sets included in the concensus tree:" << endl << endl;
		for (it2 = nodePairsInTree.begin(); it2 != nodePairsInTree.end(); it2++) {
			//only output pairs not leaves
			if (it2->first.size() > 1) { 
				temp.clear();
				//initialize temp to all "."
				temp.resize(treeSet.size(), ".");
				
				//set the spot in temp that represents it2->first[i] to a "*" 
				for (int i = 0; i < it2->first.size(); i++) {
					//find spot 
					int index = findSpot(it2->first[i]);
					temp[index] = "*";
				}
				
				//output temp
				for (int j = 0; j < temp.size(); j++) { 
					out2 << temp[j];
				}
				out2 << '\t' << it2->second << endl;
			}
		}
		
		//output sets not included
		out2 << endl << "Sets NOT included in the concensus tree:" << endl << endl;
		for (it2 = nodePairs.begin(); it2 != nodePairs.end(); it2++) {
			temp.clear();
			//initialize temp to all "."
			temp.resize(treeSet.size(), ".");
				
			//set the spot in temp that represents it2->first[i] to a "*" 
			for (int i = 0; i < it2->first.size(); i++) {
				//find spot 
				int index = findSpot(it2->first[i]);
				temp[index] = "*";
			}
				
			//output temp
			for (int j = 0; j < temp.size(); j++) { 
				out2 << temp[j];
			}
			out2 << '\t' << it2->second << endl;
		}
		
		outputFile = getRootName(globaldata->inputFileName) + "concensus.tre";
		openOutputFile(outputFile, out);
		
		concensusTree->printForBoot(out);
		
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
int ConcensusCommand::buildConcensusTree(vector<string> nodeSet) {
	try {
		vector<string> leftChildSet;
		vector<string> rightChildSet;
		
		//if you are at a leaf
		if (nodeSet.size() == 1) {
			//return the vector index of the leaf you are at
			return concensusTree->getIndex(nodeSet[0]);
		//terminate recursion
		}else if (count == numNodes) { return 0; }
		else {
			leftChildSet = getNextAvailableSet(nodeSet);
			rightChildSet = getRestSet(nodeSet, leftChildSet);
			int left = buildConcensusTree(leftChildSet);
			int right = buildConcensusTree(rightChildSet);
			concensusTree->tree[count].setChildren(left, right);
			concensusTree->tree[count].setLabel(nodePairsInTree[nodeSet]); 
			concensusTree->tree[left].setParent(count);
			concensusTree->tree[right].setParent(count);
			count++;
			return (count-1);
		}
	
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ConcensusCommand class Function buildConcensusTree. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ConcensusCommand class function buildConcensusTree. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		
}

//**********************************************************************************************************************
void ConcensusCommand::getSets() {
	try {
		vector<string> temp;
		treeSet.clear();
		
		//for each tree add the possible pairs you find
		for (int i = 0; i < t.size(); i++) {
			
			//for each non-leaf node get descendant info.
			for (int j = numLeaves; j < numNodes; j++) {
				temp.clear();
				//go through pcounts and pull out descendants
				for (it = t[i]->tree[j].pcount.begin(); it != t[i]->tree[j].pcount.end(); it++) {
					temp.push_back(it->first);
				}
				
				//sort temp
				sort(temp.begin(), temp.end());
				
				it2 = nodePairs.find(temp);
				if (it2 != nodePairs.end()) {
					nodePairs[temp]++;
				}else{
					nodePairs[temp] = 1;
				}
			}
		}
		
		//add each leaf to terminate recursion in concensus
		//you want the leaves in there but with insignifigant sightings value so it is added last
		//for each leaf node get descendant info.
		for (int j = 0; j < numLeaves; j++) {
			
			//only need the first one since leaves have no descendants but themselves
			it = t[0]->tree[j].pcount.begin(); 
			temp.clear();  temp.push_back(it->first);
			
			//fill treeSet
			treeSet.push_back(it->first);
			
			//add leaf to list but with sighting value less then all non leaf pairs 
			nodePairs[temp] = 0;
		}

		sort(treeSet.begin(), treeSet.end());
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ConcensusCommand class Function getSets. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ConcensusCommand class function getSets. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		
}

//**********************************************************************************************************************
vector<string> ConcensusCommand::getNextAvailableSet(vector<string> bigset) {
	try {
		vector<string> largest; largest.clear();
		int largestSighting = -1;
		
		//go through the sets
		for (it2 = nodePairs.begin(); it2 != nodePairs.end(); it2++) {
			//are you a subset of bigset
			if (isSubset(bigset, it2->first) == true) {
			
				//are you the largest. if you are the same size as current largest refer to sighting
				if (it2->first.size() > largest.size()) { largest = it2->first;  largestSighting = it2->second; }
				else if (it2->first.size() == largest.size()) {
					if (it2->second > largestSighting) { largest = it2->first;  largestSighting = it2->second; }
				}
				
			}
		}
		
		//save for printing out later and for branch lengths
		nodePairsInTree[largest] = nodePairs[largest];
		
		//delete whatever set you return because it is no longer available
		nodePairs.erase(largest);
		
		return largest;
	
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ConcensusCommand class Function getNextAvailableSet. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ConcensusCommand class function getNextAvailableSet. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		

}

//**********************************************************************************************************************
vector<string> ConcensusCommand::getRestSet(vector<string> bigset, vector<string> subset) {
	try {
		vector<string> rest;
		
		for (int i = 0; i < bigset.size(); i++) {
			bool inSubset = false;
			for (int j = 0; j < subset.size(); j++) {
				if (bigset[i] == subset[j]) { inSubset = true; break; }
			}
			
			//its not in the subset so put it in the rest
			if (inSubset == false) { rest.push_back(bigset[i]); }
		}
		
		//save for printing out later and for branch lengths
		nodePairsInTree[rest] = nodePairs[rest];
		
		//delete whatever set you return because it is no longer available
		nodePairs.erase(rest);

		return rest;
	
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ConcensusCommand class Function getRestSet. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ConcensusCommand class function getRestSet. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		

}

//**********************************************************************************************************************
bool ConcensusCommand::isSubset(vector<string> bigset, vector<string> subset) {
	try {
		
		//check if each guy in suset is also in bigset
		for (int i = 0; i < subset.size(); i++) {
			bool match = false;
			for (int j = 0; j < bigset.size(); j++) {
				if (subset[i] == bigset[j]) { match = true; break; }
			}
			
			//you have a guy in subset that had no match in bigset
			if (match == false) { return false; }
		}
		
		return true;
	
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ConcensusCommand class Function isSubset. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ConcensusCommand class function isSubset. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		
}
//**********************************************************************************************************************
int ConcensusCommand::findSpot(string node) {
	try {
		int spot;
		
		//check if each guy in suset is also in bigset
		for (int i = 0; i < treeSet.size(); i++) {
			if (treeSet[i] == node) { spot = i; break; }
		}
		
		return spot;
	
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ConcensusCommand class Function findSpot. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ConcensusCommand class function findSpot. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		
}
//**********************************************************************************************************************




