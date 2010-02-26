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

ConcensusCommand::ConcensusCommand(string fileroot)  {
	try {
		globaldata = GlobalData::getInstance();
		abort = false;
		
		filename = fileroot;
		
		t = globaldata->gTree;
	
	}
	catch(exception& e) {
		m->errorOut(e, "ConcensusCommand", "ConcensusCommand");
		exit(1);
	}
}

//**********************************************************************************************************************

void ConcensusCommand::help(){
	try {
		m->mothurOut("The consensus command can only be executed after a successful read.tree command.\n");
		m->mothurOut("The consensus command has no parameters.\n");
		m->mothurOut("The consensus command should be in the following format: consensus().\n");
		m->mothurOut("The consensus command output two files: .consensus.tre and .consensuspairs.\n");
		m->mothurOut("The .consensus.tre file contains the consensus tree of the trees in your input file.\n");
		m->mothurOut("The branch lengths are the percentage of trees in your input file that had the given pair.\n");
		m->mothurOut("The .consensuspairs file contains a list of the internal nodes in your tree.  For each node, the pair that was used in the consensus tree \n");
		m->mothurOut("is reported with its percentage, as well as the other pairs that were seen for that node but not used and their percentages.\n\n");		
	}
	catch(exception& e) {
		m->errorOut(e, "ConcensusCommand", "help");
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
		
		//output sets included
		out2 << endl << "Sets included in the consensus tree:" << endl << endl;
		
		vector<string> temp;
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
					//temp[index] = it2->first[i] + "  ";
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
		m->errorOut(e, "ConcensusCommand", "execute");
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
			//finds best child pair
			leftChildSet = getNextAvailableSet(nodeSet, rightChildSet);
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
		m->errorOut(e, "ConcensusCommand", "buildConcensusTree");
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
		
		
		map< vector<string>, int> nodePairsCopy = nodePairs;
		
		
		//set initial rating on pairs to sightings + subgroup sightings
		while (nodePairsCopy.size() != 0) {
		
			vector<string> small = getSmallest(nodePairsCopy);
			
			int subgrouprate = getSubgroupRating(small);
		
			nodePairsInitialRate[small] = nodePairs[small] + subgrouprate;
			
			nodePairsCopy.erase(small);
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "ConcensusCommand", "getSets");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> ConcensusCommand::getSmallest(map< vector<string>, int> nodes) {
	try{
		vector<string> smallest = nodes.begin()->first;
		int smallsize = smallest.size();
		
		for(it2 = nodes.begin(); it2 != nodes.end(); it2++) {
			if(it2->first.size() < smallsize) { smallsize = it2->first.size();  smallest = it2->first;  }
		}
		
		return smallest;
	}
	catch(exception& e) {
		m->errorOut(e, "ConcensusCommand", "getSmallest");
		exit(1);
	}
}

//**********************************************************************************************************************
vector<string> ConcensusCommand::getNextAvailableSet(vector<string> bigset, vector<string>& rest) {
	try {
//cout << "new call " << endl << endl << endl;
		vector<string> largest; largest.clear();
		rest.clear();
		
		//if you are just 2 groups
		if (bigset.size() == 2)  {   
			rest.push_back(bigset[0]);
			largest.push_back(bigset[1]);
		}else{
			rest = bestSplit[bigset][0];
			largest = bestSplit[bigset][1];
		}
		
		
		//save for printing out later and for branch lengths
		nodePairsInTree[rest] = nodePairs[rest];
		
		//delete whatever set you return because it is no longer available
		nodePairs.erase(rest);

		//save for printing out later and for branch lengths
		nodePairsInTree[largest] = nodePairs[largest];
		
		//delete whatever set you return because it is no longer available
		nodePairs.erase(largest);
		
		return largest;
	
	}
	catch(exception& e) {
		m->errorOut(e, "ConcensusCommand", "getNextAvailableSet");
		exit(1);
	}
}

/**********************************************************************************************************************/
int ConcensusCommand::getSubgroupRating(vector<string> group) {
	try {
		map< vector<string>, int>::iterator ittemp;
		map< vector< vector<string> > , int >::iterator it3;
		int rate = 0;
		
		// ***********************************************************************************//
		//1. this function must be called passing it littlest sets to biggest 
		//		since it the rating is made from your sighting plus you best splits rating
		//2. it saves the top pair to use later
		// ***********************************************************************************//

		
		if (group.size() < 3) {  return rate;  }
		
		map< vector<string>, int> possiblePairing;  //this is all the subsets of group
		
		//go through the sets
		for (it2 = nodePairs.begin(); it2 != nodePairs.end(); it2++) {
			//are you a subset of bigset, then save in possiblePairings
			if (isSubset(group, it2->first) == true) {  possiblePairing[it2->first] = it2->second;  }
		}		
	
		map< vector< vector<string> > , int > rating;
		
		while (possiblePairing.size() != 0) {
		
			it2 = possiblePairing.begin(); 
			vector<string> temprest = getRestSet(group, it2->first);
			
			//is the rest a set available in possiblePairings
			ittemp = possiblePairing.find(temprest);
			if (ittemp != possiblePairing.end()) {  //if the rest is in the possible pairings then add this pair to rating map
				vector< vector<string> > temprate;
				temprate.push_back(it2->first);  temprate.push_back(temprest);
				
				rating[temprate] = (nodePairsInitialRate[it2->first] + nodePairsInitialRate[temprest]);
				
				//erase so you dont add 1,2 and 2,1.
				possiblePairing.erase(temprest);
			}
			
			possiblePairing.erase(it2);
		}


		it3 = rating.begin();
		rate = it3->second;
		vector< vector<string> > topPair = it3->first;
		
		//choose the split with the best rating
		for (it3 = rating.begin(); it3 != rating.end(); it3++) {
			
			if (it3->second > rate) {  rate = it3->second;  topPair = it3->first;  }
		}
		
		bestSplit[group] = topPair;
		
		return rate;
	}
	catch(exception& e) {
		m->errorOut(e, "ConcensusCommand", "getSubgroupRating");
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

		return rest;
	
	}
	catch(exception& e) {
		m->errorOut(e, "ConcensusCommand", "getRestSet");
		exit(1);
	}
}

//**********************************************************************************************************************
bool ConcensusCommand::isSubset(vector<string> bigset, vector<string> subset) {
	try {
		
	
		if (subset.size() > bigset.size()) { return false;  }
		
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
		m->errorOut(e, "ConcensusCommand", "isSubset");
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
		m->errorOut(e, "ConcensusCommand", "findSpot");
		exit(1);
	}
}
//**********************************************************************************************************************




