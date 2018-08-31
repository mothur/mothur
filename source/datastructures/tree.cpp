/*
 *  tree.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 1/22/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "tree.h"

/*****************************************************************/
Tree::Tree(int num, CountTable* t, vector<string>& T) : ct(t) {
	try {
		m = MothurOut::getInstance();
		
		numLeaves = num;  
		numNodes = 2*numLeaves - 1;
        
		tree.resize(numNodes);
        Treenames = T;
	}
	catch(exception& e) {
		m->errorOut(e, "Tree", "Tree - numNodes");
		exit(1);
	}
}
/*****************************************************************/
Tree::Tree(CountTable* t, vector<string>& Tnames) : ct(t) {
	try {
		m = MothurOut::getInstance();
        
        if (Tnames.size() == 0) {  m->mothurOut("[ERROR]: no valid treenames.\n"); m->setControl_pressed(true);   }
        Treenames = Tnames;
        
		numLeaves = Treenames.size();
		numNodes = 2*numLeaves - 1;
		
		tree.resize(numNodes);
			
		//initialize groupNodeInfo
        vector<string> namesOfGroups = ct->getNamesOfGroups();
		for (int i = 0; i < namesOfGroups.size(); i++) {  groupNodeInfo[namesOfGroups[i]].resize(0);  }
		
		//initialize tree with correct number of nodes, name and group info.
		for (int i = 0; i < numNodes; i++) {
			//initialize leaf nodes
			if (i <= (numLeaves-1)) {
				tree[i].setName(Treenames[i]);
				
				//save group info
                int maxPars = 1;
				vector<string> group;
                vector<int> counts = ct->getGroupCounts(Treenames[i]);
				for (int j = 0; j < namesOfGroups.size(); j++) {  
                    if (counts[j] != 0) { //you have seqs from this group
                        groupNodeInfo[namesOfGroups[j]].push_back(i);
                        group.push_back(namesOfGroups[j]);
                        tree[i].pGroups[namesOfGroups[j]] = counts[j];
                        tree[i].pcount[namesOfGroups[j]] = counts[j];
                        //keep highest group
						if(counts[j] > maxPars){ maxPars = counts[j]; }
                    }  
                }
				tree[i].setGroup(group);
				setIndex(Treenames[i], i);
                
                if (maxPars > 1) { //then we have some more dominant groups
					//erase all the groups that are less than maxPars because you found a more dominant group.
					for(it=tree[i].pGroups.begin();it!=tree[i].pGroups.end();){
						if(it->second < maxPars){
							tree[i].pGroups.erase(it++);
						}else { it++; }
					}
					//set one remaining groups to 1
					for(it=tree[i].pGroups.begin();it!=tree[i].pGroups.end();it++){
						tree[i].pGroups[it->first] = 1;
					}
				}//end if
                
			//intialize non leaf nodes
			}else if (i > (numLeaves-1)) {
				tree[i].setName("");
				vector<string> tempGroups;
				tree[i].setGroup(tempGroups);
			}
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "Tree", "Tree");
		exit(1);
	}
}
/*****************************************************************/
Tree::Tree(CountTable* t, vector< vector<double> >& sims, vector<string>& Tnames) : ct(t) {
	try {
		m = MothurOut::getInstance();
		
		if (Tnames.size() == 0) {  m->mothurOut("[ERROR]: no valid treenames.\n"); m->setControl_pressed(true);   }
        Treenames = Tnames;
        
		numLeaves = Treenames.size();
		numNodes = 2*numLeaves - 1;
		
		tree.resize(numNodes);
        
		//initialize groupNodeInfo
        vector<string> namesOfGroups = ct->getNamesOfGroups();
		for (int i = 0; i < namesOfGroups.size(); i++) {  groupNodeInfo[namesOfGroups[i]].resize(0);  }
		
		//initialize tree with correct number of nodes, name and group info.
		for (int i = 0; i < numNodes; i++) {
			//initialize leaf nodes
			if (i <= (numLeaves-1)) {
				tree[i].setName(Treenames[i]);
				
				//save group info
                int maxPars = 1;
				vector<string> group;
                vector<int> counts = ct->getGroupCounts(Treenames[i]);
				for (int j = 0; j < namesOfGroups.size(); j++) {  
                    if (counts[j] != 0) { //you have seqs from this group
                        groupNodeInfo[namesOfGroups[j]].push_back(i);
                        group.push_back(namesOfGroups[j]);
                        tree[i].pGroups[namesOfGroups[j]] = counts[j];
                        tree[i].pcount[namesOfGroups[j]] = counts[j];
                        //keep highest group
						if(counts[j] > maxPars){ maxPars = counts[j]; }
                    }  
                }
				tree[i].setGroup(group);
				setIndex(Treenames[i], i);
                
                if (maxPars > 1) { //then we have some more dominant groups
					//erase all the groups that are less than maxPars because you found a more dominant group.
					for(it=tree[i].pGroups.begin();it!=tree[i].pGroups.end();){
						if(it->second < maxPars){
							tree[i].pGroups.erase(it++);
						}else { it++; }
					}
					//set one remaining groups to 1
					for(it=tree[i].pGroups.begin();it!=tree[i].pGroups.end();it++){
						tree[i].pGroups[it->first] = 1;
					}
				}//end if
                
                //intialize non leaf nodes
			}else if (i > (numLeaves-1)) {
				tree[i].setName("");
				vector<string> tempGroups;
				tree[i].setGroup(tempGroups);
			}
		}

        
        //build tree from matrix
        //initialize indexes
        map<int, int> thisIndexes;  //maps row in simMatrix to vector index in the tree
        for (int g = 0; g < numLeaves; g++) {	thisIndexes[g] = g;	}
		
		//do merges and create tree structure by setting parents and children
		//there are numGroups - 1 merges to do
		for (int i = 0; i < (numLeaves - 1); i++) {
			float largest = -1000.0;
			
			if (m->getControl_pressed()) { break; }
			
			int row, column;
			//find largest value in sims matrix by searching lower triangle
			for (int j = 1; j < sims.size(); j++) {
				for (int k = 0; k < j; k++) {
					if (sims[j][k] > largest) {  largest = sims[j][k]; row = j; column = k;  }
				}
			}
            
			//set non-leaf node info and update leaves to know their parents
			//non-leaf
			tree[numLeaves + i].setChildren(thisIndexes[row], thisIndexes[column]);
			
			//parents
			tree[thisIndexes[row]].setParent(numLeaves + i);
			tree[thisIndexes[column]].setParent(numLeaves + i);
			
			//blength = distance / 2;
			float blength = ((1.0 - largest) / 2);
			
			//branchlengths
			tree[thisIndexes[row]].setBranchLength(blength - tree[thisIndexes[row]].getLengthToLeaves());
			tree[thisIndexes[column]].setBranchLength(blength - tree[thisIndexes[column]].getLengthToLeaves());
			
			//set your length to leaves to your childs length plus branchlength
			tree[numLeaves + i].setLengthToLeaves(tree[thisIndexes[row]].getLengthToLeaves() + tree[thisIndexes[row]].getBranchLength());
			
			
			//update index 
			thisIndexes[row] = numLeaves+i;
			thisIndexes[column] = numLeaves+i;
			
			//remove highest value that caused the merge.
			sims[row][column] = -1000.0;
			sims[column][row] = -1000.0;
			
			//merge values in simsMatrix
			for (int n = 0; n < sims.size(); n++)	{
				//row becomes merge of 2 groups
				sims[row][n] = (sims[row][n] + sims[column][n]) / 2;
				sims[n][row] = sims[row][n];
				//delete column
				sims[column][n] = -1000.0;
				sims[n][column] = -1000.0;
			}
		}
		
		//adjust tree to make sure root to tip length is .5
		int root = findRoot();
		tree[root].setBranchLength((0.5 - tree[root].getLengthToLeaves()));
        
    }
	catch(exception& e) {
		m->errorOut(e, "Tree", "Tree");
		exit(1);
	}
}
/*****************************************************************/
Tree::~Tree() { }
/*****************************************************************/
int Tree::getIndex(string searchName) {
	try {
        map<string, int>::iterator itIndex = indexes.find(searchName);
        if (itIndex != indexes.end()) {
            return itIndex->second;
        }
		return -1;
	}
	catch(exception& e) {
		m->errorOut(e, "Tree", "getIndex");
		exit(1);
	}
}
/*****************************************************************/

void Tree::setIndex(string searchName, int index) {
	try {
		map<string, int>::iterator itIndex = indexes.find(searchName);
        if (itIndex == indexes.end()) {
            indexes[searchName] = index;
        }
	}
	catch(exception& e) {
		m->errorOut(e, "Tree", "setIndex");
		exit(1);
	}
}
/*****************************************************************/
int Tree::assembleTree() {
	try {
        //initialize groupNodeInfo
        for (int i = 0; i < (ct->getNamesOfGroups()).size(); i++) { groupNodeInfo[(ct->getNamesOfGroups())[i]].resize(0); }
        
		//build the pGroups in non leaf nodes to be used in the parsimony calcs.
		for (int i = numLeaves; i < numNodes; i++) {
			if (m->getControl_pressed()) { return 1; }

			tree[i].pGroups = (mergeGroups(i));
			tree[i].pcount = (mergeGcounts(i));
		}
        
        for(int i = 0; i < numLeaves; i++){ for (int k = 0; k < (tree[i].getGroup()).size(); k++) {  groupNodeInfo[(tree[i].getGroup())[k]].push_back(i); } }
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "Tree", "assembleTree");
		exit(1);
	}
}
/*****************************************************************/
//assumes leaf node names are in groups and no names file - used by indicator command
void Tree::getSubTree(Tree* Ctree, vector<string> Groups) {
	try {
        
        //copy Tree since we are going to destroy it
        vector<string> T = Ctree->getTreeNames();
        Tree* copy = new Tree(ct, T);
        copy->getCopy(Ctree);
        copy->assembleTree();
        
		//we want to select some of the leaf nodes to create the output tree
		//go through the input Tree starting at parents of leaves
        //initialize groupNodeInfo
        vector<string> namesOfGroups = ct->getNamesOfGroups();
		for (int i = 0; i < namesOfGroups.size(); i++) {  groupNodeInfo[namesOfGroups[i]].resize(0);  }
		
		//initialize tree with correct number of nodes, name and group info.
		for (int i = 0; i < numNodes; i++) {
			//initialize leaf nodes
			if (i <= (numLeaves-1)) {
				tree[i].setName(Groups[i]);
				
				//save group info
                int maxPars = 1;
				vector<string> group;
                vector<int> counts = ct->getGroupCounts(Groups[i]);
				for (int j = 0; j < namesOfGroups.size(); j++) {  
                    if (counts[j] != 0) { //you have seqs from this group
                        groupNodeInfo[namesOfGroups[j]].push_back(i);
                        group.push_back(namesOfGroups[j]);
                        tree[i].pGroups[namesOfGroups[j]] = counts[j];
                        tree[i].pcount[namesOfGroups[j]] = counts[j];
                        //keep highest group
						if(counts[j] > maxPars){ maxPars = counts[j]; }
                    }  
                }
				tree[i].setGroup(group);
				setIndex(Groups[i], i);
                
                if (maxPars > 1) { //then we have some more dominant groups
					//erase all the groups that are less than maxPars because you found a more dominant group.
					for(it=tree[i].pGroups.begin();it!=tree[i].pGroups.end();){
						if(it->second < maxPars){
							tree[i].pGroups.erase(it++);
						}else { it++; }
					}
					//set one remaining groups to 1
					for(it=tree[i].pGroups.begin();it!=tree[i].pGroups.end();it++){
						tree[i].pGroups[it->first] = 1;
					}
				}//end if
                
                //intialize non leaf nodes
			}else if (i > (numLeaves-1)) {
				tree[i].setName("");
				vector<string> tempGroups;
				tree[i].setGroup(tempGroups);
			}
		}
        Utils util;
		set<int> removedLeaves;
		for (int i = 0; i < copy->getNumLeaves(); i++) {
			
			if (removedLeaves.count(i) == 0) {
			
			//am I in the group
			int parent = copy->tree[i].getParent();
			
			if (parent != -1) {
				
				if (util.inUsersGroups(copy->tree[i].getName(), Groups)) {
					//find my siblings name
					int parentRC = copy->tree[parent].getRChild();
					int parentLC = copy->tree[parent].getLChild();
					
					//if I am the right child, then my sib is the left child
					int sibIndex = parentRC;
					if (parentRC == i) { sibIndex = parentLC; }
					
					string sibsName = copy->tree[sibIndex].getName();
					
					//if yes, is my sibling
					if ((util.inUsersGroups(sibsName, Groups)) || (sibsName == "")) {
						//we both are okay no trimming required
					}else{
						//i am, my sib is not, so remove sib by setting my parent to my grandparent
						int grandparent = copy->tree[parent].getParent();
						int grandparentLC = copy->tree[grandparent].getLChild();
						int grandparentRC = copy->tree[grandparent].getRChild();
						
						//whichever of my granparents children was my parent now equals me
						if (grandparentLC == parent) { grandparentLC = i; }
						else { grandparentRC = i; }
						
						copy->tree[i].setParent(grandparent);
						copy->tree[i].setBranchLength((copy->tree[i].getBranchLength()+copy->tree[parent].getBranchLength()));
						if (grandparent != -1) {
							copy->tree[grandparent].setChildren(grandparentLC, grandparentRC);
						}
						removedLeaves.insert(sibIndex);
					}
				}else{
					//find my siblings name
					int parentRC = copy->tree[parent].getRChild();
					int parentLC = copy->tree[parent].getLChild();
					
					//if I am the right child, then my sib is the left child
					int sibIndex = parentRC;
					if (parentRC == i) { sibIndex = parentLC; }
					
					string sibsName = copy->tree[sibIndex].getName();
					
					//if no is my sibling
					if ((util.inUsersGroups(sibsName, Groups)) || (sibsName == "")) {
						//i am not, but my sib is
						int grandparent = copy->tree[parent].getParent();
						int grandparentLC = copy->tree[grandparent].getLChild();
						int grandparentRC = copy->tree[grandparent].getRChild();
						
						//whichever of my granparents children was my parent now equals my sib
						if (grandparentLC == parent) { grandparentLC = sibIndex; }
						else { grandparentRC = sibIndex; }
						
						copy->tree[sibIndex].setParent(grandparent);
						copy->tree[sibIndex].setBranchLength((copy->tree[sibIndex].getBranchLength()+copy->tree[parent].getBranchLength()));
						if (grandparent != -1) {
							copy->tree[grandparent].setChildren(grandparentLC, grandparentRC);
						}
						removedLeaves.insert(i);
					}else{
						//neither of us are, so we want to eliminate ourselves and our parent
						//so set our parents sib to our great-grandparent
						int parent = copy->tree[i].getParent();
						int grandparent = copy->tree[parent].getParent();
						int parentsSibIndex;
						if (grandparent != -1) {
							int greatgrandparent = copy->tree[grandparent].getParent();
							int greatgrandparentLC, greatgrandparentRC;
							if (greatgrandparent != -1) {
								greatgrandparentLC = copy->tree[greatgrandparent].getLChild();
								greatgrandparentRC = copy->tree[greatgrandparent].getRChild();
							}
							
							int grandparentLC = copy->tree[grandparent].getLChild();
							int grandparentRC = copy->tree[grandparent].getRChild();
							
							parentsSibIndex = grandparentLC;
							if (grandparentLC == parent) { parentsSibIndex = grandparentRC; }

							//whichever of my greatgrandparents children was my grandparent
							if (greatgrandparentLC == grandparent) { greatgrandparentLC = parentsSibIndex; }
							else { greatgrandparentRC = parentsSibIndex; }
							
							copy->tree[parentsSibIndex].setParent(greatgrandparent);
							copy->tree[parentsSibIndex].setBranchLength((copy->tree[parentsSibIndex].getBranchLength()+copy->tree[grandparent].getBranchLength()));
							if (greatgrandparent != -1) {
								copy->tree[greatgrandparent].setChildren(greatgrandparentLC, greatgrandparentRC);
							}
						}else{
							copy->tree[parent].setParent(-1);
							
						}
						removedLeaves.insert(sibIndex);
						removedLeaves.insert(i);
					}
				}
			}
			}
		}
		
		int root = 0;
		for (int i = 0; i < copy->getNumNodes(); i++) {
			//you found the root
			if (copy->tree[i].getParent() == -1) { root = i; break; }
		}
        
		int nextSpot = numLeaves;
		populateNewTree(copy->tree, root, nextSpot);
        
        delete copy;
	}
	catch(exception& e) {
		m->errorOut(e, "Tree", "getSubTree");
		exit(1);
	}
}
/*****************************************************************/
int Tree::populateNewTree(vector<Node>& oldtree, int node, int& index) {
	try {
		
		if (oldtree[node].getLChild() != -1) {
			int rc = populateNewTree(oldtree, oldtree[node].getLChild(), index);
			int lc = populateNewTree(oldtree, oldtree[node].getRChild(), index);

			tree[index].setChildren(lc, rc);
			tree[rc].setParent(index);
			tree[lc].setParent(index);
			
			tree[index].setBranchLength(oldtree[node].getBranchLength());
			tree[rc].setBranchLength(oldtree[oldtree[node].getLChild()].getBranchLength());
			tree[lc].setBranchLength(oldtree[oldtree[node].getRChild()].getBranchLength());
			
			return (index++);
		}else { //you are a leaf
			int indexInNewTree = getIndex(oldtree[node].getName());
			return indexInNewTree;
		}
	}
	catch(exception& e) {
		m->errorOut(e, "Tree", "populateNewTree");
		exit(1);
	}
}
/*****************************************************************/
void Tree::getCopy(Tree* copy, bool subsample) {
	try {
        
		//for each node in the tree copy its info
		for (int i = 0; i < numNodes; i++) {
			//copy branch length
			tree[i].setBranchLength(copy->tree[i].getBranchLength());
            
			//copy parent
			tree[i].setParent(copy->tree[i].getParent());
            
			//copy children
			tree[i].setChildren(copy->tree[i].getLChild(), copy->tree[i].getRChild());
        }
        
        //build the pGroups in non leaf nodes to be used in the parsimony calcs.
		for (int i = numLeaves; i < numNodes; i++) {
			if (m->getControl_pressed()) { break; }
            
			tree[i].pGroups = (mergeGroups(i));
			tree[i].pcount = (mergeGcounts(i));
		}
	}
	catch(exception& e) {
		m->errorOut(e, "Tree", "getCopy");
		exit(1);
	}
}
/*****************************************************************/
void Tree::getCopy(Tree* copy) {
	try {
	
		//for each node in the tree copy its info
		for (int i = 0; i < numNodes; i++) {
			//copy name
			tree[i].setName(copy->tree[i].getName());
		
			//copy group
			tree[i].setGroup(copy->tree[i].getGroup());
			
			//copy branch length
			tree[i].setBranchLength(copy->tree[i].getBranchLength());
		
			//copy parent
			tree[i].setParent(copy->tree[i].getParent());
		
			//copy children
			tree[i].setChildren(copy->tree[i].getLChild(), copy->tree[i].getRChild());
		
			//copy index in node and tmap
            setIndex(copy->tree[i].getName(), getIndex(copy->tree[i].getName()));
			tree[i].setIndex(copy->tree[i].getIndex());
			
			//copy pGroups
			tree[i].pGroups = copy->tree[i].pGroups;
		
			//copy pcount
			tree[i].pcount = copy->tree[i].pcount;
		}
		
		groupNodeInfo = copy->groupNodeInfo;
		
	}
	catch(exception& e) {
		m->errorOut(e, "Tree", "getCopy");
		exit(1);
	}
}
/*****************************************************************/
//returns a map with a groupname and the number of times that group was seen in the children
//for instance if your children are white and black then it would return a map with 2 entries
// p[white] = 1 and p[black] = 1.  Now go up a level and merge that with a node who has p[white] = 1
//and you get p[white] = 2, p[black] = 1, but you erase the p[black] because you have a p value higher than 1.

map<string, int> Tree::mergeGroups(int i) {
	try {
		int lc = tree[i].getLChild();
		int rc = tree[i].getRChild();

		//set parsimony groups to left child
		map<string,int> parsimony = tree[lc].pGroups;
		
		int maxPars = 1;

		//look at right child groups and update maxPars if right child has something higher for that group.
		for(it=tree[rc].pGroups.begin();it!=tree[rc].pGroups.end();it++){
			it2 = parsimony.find(it->first);
			if (it2 != parsimony.end()) { parsimony[it->first]++;  }
			else { parsimony[it->first] = 1; }
			
			if(parsimony[it->first] > maxPars){ maxPars = parsimony[it->first]; }
		}
	
		// this is true if right child had a greater parsimony for a certain group
		if(maxPars > 1){
			//erase all the groups that are only 1 because you found something with 2.
			for(it=parsimony.begin();it!=parsimony.end();){
				if(it->second == 1){
					parsimony.erase(it++);
				}else { it++; }
			}
			//set one remaining groups to 1
			//so with our above example p[white] = 2 would be left and it would become p[white] = 1
			for(it=parsimony.begin();it!=parsimony.end();it++){ parsimony[it->first] = 1; }
		}
	
		return parsimony;
	}
	catch(exception& e) {
		m->errorOut(e, "Tree", "mergeGroups");
		exit(1);
	}
}
/*****************************************************************/
//returns a map with a groupname and the number of times that group was seen in the children
//for instance if your children are white and black then it would return a map with 2 entries
// p[white] = 1 and p[black] = 1.  Now go up a level and merge that with a node who has p[white] = 1
//and you get p[white] = 2, p[black] = 1, but you erase the p[black] because you have a p value higher than 1.

map<string, int> Tree::mergeUserGroups(int i, vector<string> g) {
	try {
		int lc = tree[i].getLChild();
		int rc = tree[i].getRChild();
		
        Utils util;
		//loop through nodes groups removing the ones the user doesn't want
		for(it=tree[lc].pGroups.begin();it!=tree[lc].pGroups.end();){
				if (util.inUsersGroups(it->first, g) != true) {
					tree[lc].pGroups.erase(it++);
				}else { it++; }
		}

		//loop through nodes groups removing the ones the user doesn't want
		for(it=tree[rc].pGroups.begin();it!=tree[rc].pGroups.end();){
				if (util.inUsersGroups(it->first, g) != true) {
					tree[rc].pGroups.erase(it++);
				}else { it++; }
		}

		//set parsimony groups to left child
		map<string,int> parsimony = tree[lc].pGroups;
		
		int maxPars = 1;

		//look at right child groups and update maxPars if right child has something higher for that group.
		for(it=tree[rc].pGroups.begin();it!=tree[rc].pGroups.end();it++){
			it2 = parsimony.find(it->first);
			if (it2 != parsimony.end()) {
				parsimony[it->first]++;
			}else {
				parsimony[it->first] = 1;
			}
			
			if(parsimony[it->first] > maxPars){
				maxPars = parsimony[it->first];
			}
		}
			
		// this is true if right child had a greater parsimony for a certain group
		if(maxPars > 1){
			//erase all the groups that are only 1 because you found something with 2.
			for(it=parsimony.begin();it!=parsimony.end();){
				if(it->second == 1){
					parsimony.erase(it++);
				}else { it++; }
			}

			for(it=parsimony.begin();it!=parsimony.end();it++){
				parsimony[it->first] = 1;
			}
		}		
		
		return parsimony;
	}
	catch(exception& e) {
		m->errorOut(e, "Tree", "mergeUserGroups");
		exit(1);
	}
}
/**************************************************************************************************/

map<string,int> Tree::mergeGcounts(int position) {
	try{
		map<string,int>::iterator pos;
	
		int lc = tree[position].getLChild();
		int rc = tree[position].getRChild();
	
		map<string,int> sum = tree[lc].pcount;
    
		for(it=tree[rc].pcount.begin();it!=tree[rc].pcount.end();it++){ sum[it->first] += it->second; }
        
		return sum;
	}
	catch(exception& e) {
		m->errorOut(e, "Tree", "mergeGcounts");
		exit(1);
	}
}
/**************************************************************************************************/
int Tree::randomLabels(vector<int>& nodesToSwap) {
    try {
        if (nodesToSwap.size() < 1)  {  return 0; } //nothing to swap
        
        for(int j = 0; j < nodesToSwap.size()-1;){
            
            if (m->getControl_pressed()) { break; }
            
            int z = nodesToSwap[j];
            int i = nodesToSwap[j+1];
    
            swapLabels(z,i);
            j += 2;
        }
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "Tree", "randomLabels");
        exit(1);
    }
}

/**************************************************************************************************/
//you only want to randomize the nodes that are from a group the user wants analyzed, so
//if either of the leaf nodes you are about to switch are not in the users groups then you don't want to switch them.
int Tree::swapLabels(int first, int second) {
    try {
        if ((first > numLeaves) || (second > numLeaves)) { m->mothurOut("[ERROR]: cannot swap tree indexes.\n"); m->setControl_pressed(true); return 0; }
        
        //switches node i and node z's info.
        map<string,int> lib_hold = tree[first].pGroups;
        tree[first].pGroups = (tree[second].pGroups);
        tree[second].pGroups = (lib_hold);
        
        vector<string> zgroup = tree[first].getGroup();
        tree[first].setGroup(tree[second].getGroup());
        tree[second].setGroup(zgroup);
        
        string zname = tree[first].getName();
        tree[first].setName(tree[second].getName());
        setIndex(tree[second].getName(), first);
        tree[second].setName(zname);
        setIndex(zname, second);
        
        map<string,int> gcount_hold = tree[first].pcount;
        tree[first].pcount = (tree[second].pcount);
        tree[second].pcount = (gcount_hold);
        
        return 1;
    }
    catch(exception& e) {
        m->errorOut(e, "Tree", "swapLabels");
        exit(1);
    }
}
/*************************************************************************************************/
void Tree::assembleRandomUnifracTree(vector<int> g) {
	randomLabels(g);
	assembleTree();
}
/*************************************************************************************************/
//for now it's just random topology but may become random labels as well later that why this is such a simple function now...
void Tree::assembleRandomTree(Utils* myUtil) {
	randomTopology(myUtil);
	assembleTree();
}
/**************************************************************************************************/

void Tree::randomTopology(Utils* myUtil) {
	try {
        for(int i=0;i<numNodes;i++)         { tree[i].setParent(-1);        }
        for(int i=numLeaves;i<numNodes;i++) { tree[i].setChildren(-1, -1);  }
        
        for(int i=numLeaves;i<numNodes;i++){
            int escape =0;
            int rnd_index1, rnd_index2;
            while(escape == 0){
                rnd_index1 = myUtil->getRandomIndex(i);
                if(tree[rnd_index1].getParent() == -1){escape = 1;}
            }
            
            escape = 0;
            while(escape == 0){
                rnd_index2 = myUtil->getRandomIndex(i);
                if(rnd_index2 != rnd_index1 && tree[rnd_index2].getParent() == -1){
                    escape = 1;
                }		
            }
            
            tree[i].setChildren(rnd_index1,rnd_index2);
            tree[i].setParent(-1);
            tree[rnd_index1].setParent(i);
            tree[rnd_index2].setParent(i);
        }
	}
	catch(exception& e) {
		m->errorOut(e, "Tree", "randomTopology");
		exit(1);
	}
}
/*****************************************************************/
vector<int> Tree::getNodes(vector<string> theseGroups) {
    try {
        set<int> nodes;
        for (int i = 0; i < theseGroups.size(); i++) {
            if (m->getControl_pressed()) { break; }
            
            map<string, vector<int> >::iterator it = groupNodeInfo.find(theseGroups[i]);
            if (it != groupNodeInfo.end()) {//we have nodes for this group
                for (int j = 0; j < (it->second).size(); j++) { nodes.insert((it->second)[j]); } //removes dups
            }
        }
        
        vector<int> uniqueNodes;
        for (set<int>::iterator it = nodes.begin(); it != nodes.end(); it++) { uniqueNodes.push_back(*it); }
        
        return uniqueNodes;
    }
    catch(exception& e) {
        m->errorOut(e, "Tree", "getNodes");
        exit(1);
    }
}
/*****************************************************************/
void Tree::print(ostream& out) {
	try {
		int root = findRoot();
		printBranch(root, out, "branch");
		out << ";" << endl;
	}
	catch(exception& e) {
		m->errorOut(e, "Tree", "print");
		exit(1);
	}
}
/*****************************************************************/
void Tree::print(ostream& out, map<string, string> nameMap) {
	try {
		int root = findRoot();
		printBranch(root, out, nameMap);
		out << ";" << endl;
	}
	catch(exception& e) {
		m->errorOut(e, "Tree", "print");
		exit(1);
	}
}
/*****************************************************************/
void Tree::print(ostream& out, string mode) {
	try {
		int root = findRoot();
		printBranch(root, out, mode);
		out << ";" << endl;
	}
	catch(exception& e) {
		m->errorOut(e, "Tree", "print");
		exit(1);
	}
}
/*****************************************************************/
// This prints out the tree in Newick form.
void Tree::createNewickFile(string f) {
	try {
		int root = findRoot();
	
		filename = f;

        Utils util; util.openOutputFile(filename, out);
		
		printBranch(root, out, "branch");
		
		// you are at the end of the tree
		out << ";" << endl;
		out.close();
	}
	catch(exception& e) {
		m->errorOut(e, "Tree", "createNewickFile");
		exit(1);
	}
}

/*****************************************************************/
//This function finds the index of the root node.

int Tree::findRoot() {
	try {
		for (int i = 0; i < numNodes; i++) {
			//you found the root
			if (tree[i].getParent() == -1) { return i; }
			
		}
		return -1;
	}
	catch(exception& e) {
		m->errorOut(e, "Tree", "findRoot");
		exit(1);
	}
}
/*****************************************************************/
void Tree::printBranch(int node, ostream& out, map<string, string> names) {
try {
// you are not a leaf
		if (tree[node].getLChild() != -1) {
			out << "(";
			printBranch(tree[node].getLChild(), out, names);
			out << ",";
			printBranch(tree[node].getRChild(), out, names);
			out << ")";
			
            //if there is a branch length then print it
            if (tree[node].getBranchLength() != -1) {
                out << ":" << tree[node].getBranchLength();
            }
			
		}else { //you are a leaf
            map<string, string>::iterator itNames = names.find(tree[node].getName());
            Utils util;
            string outputString = "";
            if (itNames != names.end()) { 
                
                vector<string> dupNames;
                util.splitAtComma((itNames->second), dupNames);
                
                if (dupNames.size() == 1) {
                    outputString += tree[node].getName();
                    if (tree[node].getBranchLength() != -1) {
                        outputString += ":" + toString(tree[node].getBranchLength());
                    }
                }else {
                    outputString += "(";
                    
                    for (int u = 0; u < dupNames.size()-1; u++) {
                        outputString += dupNames[u];
                        
                        if (tree[node].getBranchLength() != -1) {
                            outputString += ":" + toString(0.0);
                        }
                        outputString += ",";
                    }
                    
                    outputString += dupNames[dupNames.size()-1];
                    if (tree[node].getBranchLength() != -1) {
                        outputString += ":" + toString(0.0);
                    }
                    
                    outputString += ")";
                    if (tree[node].getBranchLength() != -1) {
                        outputString += ":" + toString(tree[node].getBranchLength());
                    }
                }
            }else { 
                outputString = tree[node].getName();
                //if there is a branch length then print it
                if (tree[node].getBranchLength() != -1) {
                    outputString += ":" + toString(tree[node].getBranchLength());
                }
                
                m->mothurOut("[ERROR]: " + tree[node].getName() + " is not in your namefile, please correct."); m->mothurOutEndLine(); 
            }
                
            out << outputString;
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "Tree", "printBranch");
		exit(1);
	}
}
/*****************************************************************/
void Tree::printBranch(int node, ostream& out, string mode) {
    try {
        
        // you are not a leaf
		if (tree[node].getLChild() != -1) {
			out << "(";
			printBranch(tree[node].getLChild(), out, mode);
			out << ",";
			printBranch(tree[node].getRChild(), out, mode);
			out << ")";
			if (mode == "branch") {
				//if there is a branch length then print it
				if (tree[node].getBranchLength() != -1) {
					out << ":" << tree[node].getBranchLength();
				}
			}else if (mode == "boot") {
				//if there is a label then print it
				if (tree[node].getLabel() != "") {
					out << tree[node].getLabel();
				}
			}else if (mode == "both") {
				if (tree[node].getLabel() != "") {
					out << tree[node].getLabel();
				}
				//if there is a branch length then print it
				if (tree[node].getBranchLength() != -1) {
					out << ":" << tree[node].getBranchLength();
				}
			}
		}else { //you are a leaf
			vector<string> leafGroup = ct->getGroups(tree[node].getName());
			
			if (mode == "branch") {
				out << leafGroup[0]; 
				//if there is a branch length then print it
				if (tree[node].getBranchLength() != -1) {
					out << ":" << tree[node].getBranchLength();
				}
			}else if (mode == "boot") {
				out << leafGroup[0]; 
				//if there is a label then print it
				if (tree[node].getLabel() != "") {
					out << tree[node].getLabel();
				}
			}else if (mode == "both") {
				out << tree[node].getName();
				if (tree[node].getLabel() != "") {
					out << tree[node].getLabel();
				}
				//if there is a branch length then print it
				if (tree[node].getBranchLength() != -1) {
					out << ":" << tree[node].getBranchLength();
				}
			}
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "Tree", "printBranch");
		exit(1);
	}
}
/*****************************************************************/
void Tree::printBranch(int node, ostream& out, string mode, vector<Node>& theseNodes) {
	try {
		
		// you are not a leaf
		if (theseNodes[node].getLChild() != -1) {
			out << "(";
			printBranch(theseNodes[node].getLChild(), out, mode);
			out << ",";
			printBranch(theseNodes[node].getRChild(), out, mode);
			out << ")";
			if (mode == "branch") {
				//if there is a branch length then print it
				if (theseNodes[node].getBranchLength() != -1) {
					out << ":" << theseNodes[node].getBranchLength();
				}
			}else if (mode == "boot") {
				//if there is a label then print it
				if (theseNodes[node].getLabel() != "") {
					out << theseNodes[node].getLabel();
				}
			}else if (mode == "both") {
				if (theseNodes[node].getLabel() != "") {
					out << theseNodes[node].getLabel();
				}
				//if there is a branch length then print it
				if (theseNodes[node].getBranchLength() != -1) {
					out << ":" << theseNodes[node].getBranchLength();
				}
			}
		}else { //you are a leaf
			vector<string> leafGroup = ct->getGroups(theseNodes[node].getName());
			
			if (mode == "branch") {
				out << leafGroup[0]; 
				//if there is a branch length then print it
				if (theseNodes[node].getBranchLength() != -1) {
					out << ":" << theseNodes[node].getBranchLength();
				}
			}else if (mode == "boot") {
				out << leafGroup[0]; 
				//if there is a label then print it
				if (theseNodes[node].getLabel() != "") {
					out << theseNodes[node].getLabel();
				}
			}else if (mode == "both") {
				out << theseNodes[node].getName();
				if (theseNodes[node].getLabel() != "") {
					out << theseNodes[node].getLabel();
				}
				//if there is a branch length then print it
				if (theseNodes[node].getBranchLength() != -1) {
					out << ":" << theseNodes[node].getBranchLength();
				}
			}
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "Tree", "printBranch");
		exit(1);
	}
}
/*****************************************************************/
void Tree::printTree() {
	for(int i=0;i<numNodes;i++){
		cout << i << '\t';
		tree[i].printNode();
	}
}
/*******************************************************/

