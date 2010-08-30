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
Tree::Tree(string g) {
	try {
		globaldata = GlobalData::getInstance();
		m = MothurOut::getInstance();
		
		parseTreeFile();  globaldata->runParse = false;  
	}
	catch(exception& e) {
		m->errorOut(e, "Tree", "Tree - just parse");
		exit(1);
	}
}
/*****************************************************************/
Tree::Tree() {
	try {
		globaldata = GlobalData::getInstance();
		m = MothurOut::getInstance();
		
		if (globaldata->runParse == true) {  parseTreeFile();  globaldata->runParse = false;  }
//for(int i = 0; i < 	globaldata->Treenames.size(); i++) { cout << i << '\t' << globaldata->Treenames[i] << endl;  }	
		numLeaves = globaldata->Treenames.size();
		numNodes = 2*numLeaves - 1;
		
		tree.resize(numNodes);
		
		//initialize groupNodeInfo
		for (int i = 0; i < globaldata->gTreemap->namesOfGroups.size(); i++) {
			groupNodeInfo[globaldata->gTreemap->namesOfGroups[i]].resize(0);
		}

		//initialize tree with correct number of nodes, name and group info.
		for (int i = 0; i < numNodes; i++) {
			//initialize leaf nodes
			if (i <= (numLeaves-1)) {
				tree[i].setName(globaldata->Treenames[i]);
				
				//save group info
				string group = globaldata->gTreemap->getGroup(globaldata->Treenames[i]);
				vector<string> tempGroups; tempGroups.push_back(group);
				tree[i].setGroup(tempGroups);
				groupNodeInfo[group].push_back(i); 
				
				//set pcount and pGroup for groupname to 1.
				tree[i].pcount[group] = 1;
				tree[i].pGroups[group] = 1;
				
				//Treemap knows name, group and index to speed up search
				globaldata->gTreemap->setIndex(globaldata->Treenames[i], i);
	
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
Tree::~Tree() {}
/*****************************************************************/
void Tree::addNamesToCounts() {
	try {
		//ex. seq1	seq2,seq3,se4
		//		seq1 = pasture
		//		seq2 = forest
		//		seq4 = pasture
		//		seq3 = ocean
		
		//before this function seq1.pcount = pasture -> 1
		//after				   seq1.pcount = pasture -> 2, forest -> 1, ocean -> 1
		
		//before this function seq1.pgroups = pasture -> 1
		//after				   seq1.pgroups = pasture -> 1 since that is the dominant group

				
		//go through each leaf and update its pcounts and pgroups
		
		//float A = clock();

		for (int i = 0; i < numLeaves; i++) {

			string name = tree[i].getName();
		
			map<string, string>::iterator itNames = globaldata->names.find(name);
		
			if (itNames == globaldata->names.end()) { m->mothurOut(name + " is not in your name file, please correct."); m->mothurOutEndLine(); exit(1);  }
			else {
				vector<string> dupNames;
				m->splitAtComma(globaldata->names[name], dupNames);
				
				map<string, int>::iterator itCounts;
				int maxPars = 1;
				set<string> groupsAddedForThisNode;
				for (int j = 0; j < dupNames.size(); j++) {
					
					string group = globaldata->gTreemap->getGroup(dupNames[j]);
					
					if (dupNames[j] != name) {//you already added yourself in the constructor
				
						if (groupsAddedForThisNode.count(group) == 0)  {  groupNodeInfo[group].push_back(i);  groupsAddedForThisNode.insert(group);  } //if you have not already added this node for this group, then add it
						
						//update pcounts
						itCounts = tree[i].pcount.find(group);
						if (itCounts == tree[i].pcount.end()) { //new group, add it
							tree[i].pcount[group] = 1;
						}else {
							tree[i].pcount[group]++;
						}
							
						//update pgroups
						itCounts = tree[i].pGroups.find(group);
						if (itCounts == tree[i].pGroups.end()) { //new group, add it
							tree[i].pGroups[group] = 1;
						}else {
							tree[i].pGroups[group]++;
						}
						
						//keep highest group
						if(tree[i].pGroups[group] > maxPars){
							maxPars = tree[i].pGroups[group];
						}
					}else {	 groupsAddedForThisNode.insert(group);  } //add it so you don't add it to groupNodeInfo again
				}//end for
				
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
				
				//update groups to reflect all the groups this node represents
				vector<string> nodeGroups;
				map<string, int>::iterator itGroups;
				for (itGroups = tree[i].pcount.begin(); itGroups != tree[i].pcount.end(); itGroups++) {
					nodeGroups.push_back(itGroups->first);
				}
				tree[i].setGroup(nodeGroups);
				
			}//end else
		}//end for		
		
		//float B = clock();
		//cout << "addNamesToCounts\t" << (B - A) / CLOCKS_PER_SEC << endl;	

	}
	catch(exception& e) {
		m->errorOut(e, "Tree", "addNamesToCounts");
		exit(1);
	}
}
/*****************************************************************/
int Tree::getIndex(string searchName) {
	try {
		//Treemap knows name, group and index to speed up search
		// getIndex function will return the vector index or -1 if seq is not found.
		int index = globaldata->gTreemap->getIndex(searchName);
		return index;
		
	}
	catch(exception& e) {
		m->errorOut(e, "Tree", "getIndex");
		exit(1);
	}
}
/*****************************************************************/

void Tree::setIndex(string searchName, int index) {
	try {
		//set index in treemap
		globaldata->gTreemap->setIndex(searchName, index);
	}
	catch(exception& e) {
		m->errorOut(e, "Tree", "setIndex");
		exit(1);
	}
}
/*****************************************************************/
int Tree::assembleTree() {
	try {
		//float A = clock();

		//if user has given a names file we want to include that info in the pgroups and pcount info.
		if(globaldata->names.size() != 0) {  addNamesToCounts();  }
		
		//build the pGroups in non leaf nodes to be used in the parsimony calcs.
		for (int i = numLeaves; i < numNodes; i++) {
			if (m->control_pressed) { return 1; }

			tree[i].pGroups = (mergeGroups(i));
			tree[i].pcount = (mergeGcounts(i));
		}
		//float B = clock();
		//cout << "assembleTree\t" << (B-A) / CLOCKS_PER_SEC << endl;
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "Tree", "assembleTree");
		exit(1);
	}
}
/*****************************************************************/
int Tree::assembleTree(string n) {
	try {
		
		//build the pGroups in non leaf nodes to be used in the parsimony calcs.
		for (int i = numLeaves; i < numNodes; i++) {
			if (m->control_pressed) { return 1; }

			tree[i].pGroups = (mergeGroups(i));
			tree[i].pcount = (mergeGcounts(i));
		}
		//float B = clock();
		//cout << "assembleTree\t" << (B-A) / CLOCKS_PER_SEC << endl;
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "Tree", "assembleTree");
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
			tree[i].setIndex(copy->tree[i].getIndex());
			setIndex(copy->tree[i].getName(), getIndex(copy->tree[i].getName()));
			
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
			//set one remaining groups to 1
			//so with our above example p[white] = 2 would be left and it would become p[white] = 1
			for(it=parsimony.begin();it!=parsimony.end();it++){
				parsimony[it->first] = 1;
			}
		
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
		
		//loop through nodes groups removing the ones the user doesn't want
		for(it=tree[lc].pGroups.begin();it!=tree[lc].pGroups.end();){
				if (m->inUsersGroups(it->first, g) != true) {
					tree[lc].pGroups.erase(it++);
				}else { it++; }
		}

		//loop through nodes groups removing the ones the user doesn't want
		for(it=tree[rc].pGroups.begin();it!=tree[rc].pGroups.end();){
				if (m->inUsersGroups(it->first, g) != true) {
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
    
		for(it=tree[rc].pcount.begin();it!=tree[rc].pcount.end();it++){
			sum[it->first] += it->second;
		}
		return sum;
	}
	catch(exception& e) {
		m->errorOut(e, "Tree", "mergeGcounts");
		exit(1);
	}
}
/**************************************************************************************************/

void Tree::randomLabels(vector<string> g) {
	try {
	
		//initialize groupNodeInfo
		for (int i = 0; i < globaldata->gTreemap->namesOfGroups.size(); i++) {
			groupNodeInfo[globaldata->gTreemap->namesOfGroups[i]].resize(0);
		}
		
		for(int i = 0; i < numLeaves; i++){
			int z;
			//get random index to switch with
			z = int((float)(i+1) * (float)(rand()) / ((float)RAND_MAX+1.0));	
			
			//you only want to randomize the nodes that are from a group the user wants analyzed, so
			//if either of the leaf nodes you are about to switch are not in the users groups then you don't want to switch them.
			bool treez, treei;
		
			treez = m->inUsersGroups(tree[z].getGroup(), g);
			treei = m->inUsersGroups(tree[i].getGroup(), g);
			
			if ((treez == true) && (treei == true)) {
				//switches node i and node z's info.
				map<string,int> lib_hold = tree[z].pGroups;
				tree[z].pGroups = (tree[i].pGroups);
				tree[i].pGroups = (lib_hold);
				
				vector<string> zgroup = tree[z].getGroup();
				tree[z].setGroup(tree[i].getGroup());
				tree[i].setGroup(zgroup);
				
				string zname = tree[z].getName();
				tree[z].setName(tree[i].getName());
				tree[i].setName(zname);
				
				map<string,int> gcount_hold = tree[z].pcount;
				tree[z].pcount = (tree[i].pcount);
				tree[i].pcount = (gcount_hold);
			}
			
			for (int k = 0; k < (tree[i].getGroup()).size(); k++) {  groupNodeInfo[(tree[i].getGroup())[k]].push_back(i); }
			for (int k = 0; k < (tree[z].getGroup()).size(); k++) {  groupNodeInfo[(tree[z].getGroup())[k]].push_back(z); }
		}
	}
	catch(exception& e) {
		m->errorOut(e, "Tree", "randomLabels");
		exit(1);
	}
}
/**************************************************************************************************

void Tree::randomLabels(string groupA, string groupB) {
	try {
		int numSeqsA = globaldata->gTreemap->seqsPerGroup[groupA];
		int numSeqsB = globaldata->gTreemap->seqsPerGroup[groupB];

		vector<string> randomGroups(numSeqsA+numSeqsB, groupA);
		for(int i=numSeqsA;i<randomGroups.size();i++){
			randomGroups[i] = groupB;
		}
		random_shuffle(randomGroups.begin(), randomGroups.end());
				
		int randomCounter = 0;				
		for(int i=0;i<numLeaves;i++){
			if(tree[i].getGroup() == groupA || tree[i].getGroup() == groupB){
				tree[i].setGroup(randomGroups[randomCounter]);
				tree[i].pcount.clear();
				tree[i].pcount[randomGroups[randomCounter]] = 1;
				tree[i].pGroups.clear();
				tree[i].pGroups[randomGroups[randomCounter]] = 1;
				randomCounter++;
			}
		}
	}		
	catch(exception& e) {
		m->errorOut(e, "Tree", "randomLabels");
		exit(1);
	}
}
/**************************************************************************************************/
void Tree::randomBlengths()  {
	try {
		for(int i=numNodes-1;i>=0;i--){
			int z = int((float)(i+1) * (float)(rand()) / ((float)RAND_MAX+1.0));	
		
			float bl_hold = tree[z].getBranchLength();
			tree[z].setBranchLength(tree[i].getBranchLength());
			tree[i].setBranchLength(bl_hold);
		}
	}
	catch(exception& e) {
		m->errorOut(e, "Tree", "randomBlengths");
		exit(1);
	}
}
/*************************************************************************************************/
void Tree::assembleRandomUnifracTree(vector<string> g) {
	randomLabels(g);
	assembleTree("noNameCounts");
}
/*************************************************************************************************/
void Tree::assembleRandomUnifracTree(string groupA, string groupB) {

	vector<string> temp; temp.push_back(groupA); temp.push_back(groupB);
	randomLabels(temp);
	assembleTree("noNameCounts");
}

/*************************************************************************************************/
//for now it's just random topology but may become random labels as well later that why this is such a simple function now...
void Tree::assembleRandomTree() {
	randomTopology();
	assembleTree();
}
/**************************************************************************************************/

void Tree::randomTopology() {
	try {
		for(int i=0;i<numNodes;i++){
			tree[i].setParent(-1);
		}
		for(int i=numLeaves;i<numNodes;i++){
			tree[i].setChildren(-1, -1); 
		}
    
		for(int i=numLeaves;i<numNodes;i++){
			int escape =0;
			int rnd_index1, rnd_index2;
			while(escape == 0){
				rnd_index1 = (int)(((double)rand() / (double) RAND_MAX)*i);
				if(tree[rnd_index1].getParent() == -1){escape = 1;}
			}
		
			escape = 0;
			while(escape == 0){
				rnd_index2 = (int)(((double)rand() / (double) RAND_MAX)*i);
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
void Tree::printForBoot(ostream& out) {
	try {
		int root = findRoot();
		printBranch(root, out, "boot");
		out << ";" << endl;
	}
	catch(exception& e) {
		m->errorOut(e, "Tree", "printForBoot");
		exit(1);
	}
}

/*****************************************************************/
// This prints out the tree in Newick form.
void Tree::createNewickFile(string f) {
	try {
		int root = findRoot();
		//filename = m->getRootName(globaldata->getTreeFile()) + "newick";
		filename = f;

		m->openOutputFile(filename, out);
		
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
			//cout << "i = " << i << endl;
			//cout << "i's parent = " << tree[i].getParent() << endl;  
		}
		return -1;
	}
	catch(exception& e) {
		m->errorOut(e, "Tree", "findRoot");
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
				if (tree[node].getLabel() != -1) {
					out << tree[node].getLabel();
				}
			}
		}else { //you are a leaf
			string leafGroup = globaldata->gTreemap->getGroup(tree[node].getName());
			
			out << leafGroup; 
			if (mode == "branch") {
				//if there is a branch length then print it
				if (tree[node].getBranchLength() != -1) {
					out << ":" << tree[node].getBranchLength();
				}
			}else if (mode == "boot") {
				//if there is a label then print it
				if (tree[node].getLabel() != -1) {
					out << tree[node].getLabel();
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

/*****************************************************************/
//this code is a mess and should be rethought...-slw
void Tree::parseTreeFile() {
	
	//only takes names from the first tree and assumes that all trees use the same names.
	try {
		string filename = globaldata->getTreeFile();
		ifstream filehandle;
		m->openInputFile(filename, filehandle);
		int c, comment;
		comment = 0;
		int done = 1;
		
		//ifyou are not a nexus file 
		if((c = filehandle.peek()) != '#') {  
			while((c = filehandle.peek()) != ';') { 
				while ((c = filehandle.peek()) != ';') {
					// get past comments
					if(c == '[') {
						comment = 1;
					}
					if(c == ']'){
						comment = 0;
					}
					if((c == '(') && (comment != 1)){ break; }
					filehandle.get();
				}

				done = readTreeString(filehandle); 
				if (done == 0) { break; }
			}
		//ifyou are a nexus file
		}else if((c = filehandle.peek()) == '#') {
			string holder = "";
					
			// get past comments
			while(holder != "translate" && holder != "Translate"){	
				if(holder == "[" || holder == "[!"){
					comment = 1;
				}
				if(holder == "]"){
					comment = 0;
				}
				filehandle >> holder; 

				//if there is no translate then you must read tree string otherwise use translate to get names
				if((holder == "tree") && (comment != 1)){	
					//pass over the "tree rep.6878900 = "
					while (((c = filehandle.get()) != '(') && ((c = filehandle.peek()) != EOF)) {;}

					if(c == EOF) { break; }
					filehandle.putback(c);  //put back first ( of tree.
					done = readTreeString(filehandle);
	
					break;
				}
			
				if (done == 0) { break;  }
			}
			
			//use nexus translation rather than parsing tree to save time
			if((holder == "translate") || (holder == "Translate")) {

				string number, name, h;
				h = ""; // so it enters the loop the first time
				while((h != ";") && (number != ";")) { 
					filehandle >> number;
					filehandle >> name;
	
					//c = , until done with translation then c = ;
					h = name.substr(name.length()-1, name.length()); 
					name.erase(name.end()-1);  //erase the comma
					globaldata->Treenames.push_back(number);
				}
				if(number == ";") { globaldata->Treenames.pop_back(); }  //in case ';' from translation is on next line instead of next to last name
			}
		}
		filehandle.close();
		
		//for (int i = 0; i < globaldata->Treenames.size(); i++) {
//cout << globaldata->Treenames[i] << endl; }
//cout << globaldata->Treenames.size() << endl;
	}
	catch(exception& e) {
		m->errorOut(e, "Tree", "parseTreeFile");
		exit(1);
	}
}
/*******************************************************/

/*******************************************************/
int Tree::readTreeString(ifstream& filehandle)	{
	try {
		int c;
		string name;  //, k
		
		while((c = filehandle.peek()) != ';') { 
//k = c;
//cout << " at beginning of while " <<  k << endl;			
			if(c == ')')  {    
				//to pass over labels in trees
				c=filehandle.get();
				while((c!=',') && (c != -1) && (c!= ':') && (c!=';')){ c=filehandle.get(); }
				filehandle.putback(c);
			}
			if(c == ';') { return 0; }
			if(c == -1) { return 0; }
			//if you are a name
			if((c != '(') && (c != ')') && (c != ',') && (c != ':') && (c != '\n') && (c != '\t') && (c != 32)) { //32 is space
				name = "";
				c = filehandle.get();
			//k = c;
//cout << k << endl;
				while ((c != '(') && (c != ')') && (c != ',') && (c != ':') && (c != '\n') && (c != 32) && (c != '\t')) {			
					name += c;
					c = filehandle.get();
			//k = c;
//cout << " in name while " << k << endl;
				}
				
//cout << "name = " << name << endl;
				globaldata->Treenames.push_back(name);
				filehandle.putback(c);
//k = c;
//cout << " after putback" <<  k << endl;
			} 
			
			if(c  == ':') { //read until you reach the end of the branch length
				while ((c != '(') && (c != ')') && (c != ',') && (c != ';') && (c != '\n') && (c != '\t') && (c != 32)) {
					c = filehandle.get();
	//k = c;
	//cout << " in branch while " << k << endl;
				}
				filehandle.putback(c);
			}
		
			c = filehandle.get();
//k = c;
	//cout << " here after get " << k << endl;
			if(c == ';') { return 0; }
			if(c == ')') { filehandle.putback(c); }
	//k = c;
//cout << k << endl;

		}
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "Tree", "readTreeString");
		exit(1);
	}
}	

/*******************************************************/

/*******************************************************/

