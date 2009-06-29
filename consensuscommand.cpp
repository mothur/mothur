/*
 *  consensuscommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 4/29/09.
 *  Copyright 2009 Schloss Lab UMASS AMherst. All rights reserved.
 *
 */

#include "consensuscommand.h"

//**********************************************************************************************************************

ConcensusCommand::ConcensusCommand(string fileroot){
	try {
		globaldata = GlobalData::getInstance();
		abort = false;
		
		filename = fileroot;
		//allow user to run help
		//if(option == "help") { help(); abort = true; }
		
		//else {
			//if (option != "") { mothurOut("There are no valid parameters for the consensus command."); mothurOutEndLine(); abort = true; }
			
		//	//no trees were read
		//	if (globaldata->gTree.size() == 0) {  mothurOut("You must execute the read.tree command, before you may use the consensus command."); mothurOutEndLine(); abort = true;  }
			//else { 
			 t = globaldata->gTree;
			 //	}
		//}
	}
	catch(exception& e) {
		errorOut(e, "ConcensusCommand", "ConcensusCommand");
		exit(1);
	}
}

//**********************************************************************************************************************

void ConcensusCommand::help(){
	try {
		mothurOut("The consensus command can only be executed after a successful read.tree command.\n");
		mothurOut("The consensus command has no parameters.\n");
		mothurOut("The consensus command should be in the following format: consensus().\n");
		mothurOut("The consensus command output two files: .consensus.tre and .consensuspairs.\n");
		mothurOut("The .consensus.tre file contains the consensus tree of the trees in your input file.\n");
		mothurOut("The branch lengths are the percentage of trees in your input file that had the given pair.\n");
		mothurOut("The .consensuspairs file contains a list of the internal nodes in your tree.  For each node, the pair that was used in the consensus tree \n");
		mothurOut("is reported with its percentage, as well as the other pairs that were seen for that node but not used and their percentages.\n\n");		
	}
	catch(exception& e) {
		errorOut(e, "ConcensusCommand", "help");
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
		
		//open file for pairing not included in the tree
		notIncluded = filename + ".cons.pairs";
		openOutputFile(notIncluded, out2);
		
		consensusTree = new Tree();
		
		it2 = nodePairs.find(treeSet);
		
		nodePairsInTree[treeSet] = it2->second; 
		
		//erase treeset because you are adding it
		nodePairs.erase(treeSet);
		
		//set count to numLeaves;
		count = numLeaves;
		
		buildConcensusTree(treeSet);
		
		consensusTree->assembleTree();
		
		//output species in order
		out2 << "Species in Order: " << endl << endl;
		for (int i = 0; i < treeSet.size(); i++) {  out2 << i+1 << ".  " << treeSet[i] << endl; }
		
		vector<string> temp; 
		//output sets included
		out2 << endl << "Sets included in the consensus tree:" << endl << endl;
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
		out2 << endl << "Sets NOT included in the consensus tree:" << endl << endl;
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
		
		outputFile = filename + ".cons.tre";
		openOutputFile(outputFile, out);
		
		consensusTree->printForBoot(out);
		
		out.close(); out2.close();
		
		delete consensusTree; 
		
		return 0;
	}
	catch(exception& e) {
		errorOut(e, "ConcensusCommand", "execute");
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
			return consensusTree->getIndex(nodeSet[0]);
		//terminate recursion
		}else if (count == numNodes) { return 0; }
		else {
			leftChildSet = getNextAvailableSet(nodeSet);
			rightChildSet = getRestSet(nodeSet, leftChildSet);
			int left = buildConcensusTree(leftChildSet);
			int right = buildConcensusTree(rightChildSet);
			consensusTree->tree[count].setChildren(left, right);
			consensusTree->tree[count].setLabel(nodePairsInTree[nodeSet]); 
			consensusTree->tree[left].setParent(count);
			consensusTree->tree[right].setParent(count);
			count++;
			return (count-1);
		}
	
	}
	catch(exception& e) {
		errorOut(e, "ConcensusCommand", "buildConcensusTree");
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
		
		//add each leaf to terminate recursion in consensus
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
		errorOut(e, "ConcensusCommand", "getSets");
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
		errorOut(e, "ConcensusCommand", "getNextAvailableSet");
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
		errorOut(e, "ConcensusCommand", "getRestSet");
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
		errorOut(e, "ConcensusCommand", "isSubset");
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
		errorOut(e, "ConcensusCommand", "findSpot");
		exit(1);
	}
}
//**********************************************************************************************************************




