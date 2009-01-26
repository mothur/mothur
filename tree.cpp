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
Tree::Tree() {
	try {
		globaldata = GlobalData::getInstance();

		numLeaves = globaldata->gTreemap->getNumSeqs();
		numNodes = 2*numLeaves - 1;
		
		tree.resize(numNodes);

		//initialize tree with correct number of nodes, name and group info.
		for (int i = 0; i < numNodes; i++) {

			//initialize leaf nodes
			if (i <= (numLeaves-1)) {
				tree[i].setName(globaldata->gTreemap->namesOfSeqs[i]);
				tree[i].setGroup(globaldata->gTreemap->getGroup(globaldata->gTreemap->namesOfSeqs[i]));
				//the node knows its index
				tree[i].setIndex(i);
				//Treemap knows name, group and index to speed up search
				globaldata->gTreemap->setIndex(globaldata->gTreemap->namesOfSeqs[i], i);
			//intialize non leaf nodes
			}else if (i > (numLeaves-1)) {
				tree[i].setName("");
				tree[i].setGroup("");
				//the node knows its index
				tree[i].setIndex(i);
			}
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the Tree class Function Tree. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the Tree class function Tree. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
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
		cout << "Standard Error: " << e.what() << " has occurred in the Tree class Function getIndex. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the Tree class function getIndex. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
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
		cout << "Standard Error: " << e.what() << " has occurred in the Tree class Function setIndex. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the Tree class function setIndex. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		
}


/*****************************************************************/
// This prints out the tree in Newick form.
void Tree::createNewickFile() {
	try {
		int root = findRoot();
		filename = getRootName(globaldata->getTreeFile()) + "newick";
		openOutputFile(filename, out);
		
		printBranch(root);
		
		// you are at the end of the tree
		out << ";" << endl;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the Tree class Function createNewickFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the Tree class function createNewickFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
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
		cout << "Standard Error: " << e.what() << " has occurred in the Tree class Function findRoot. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the Tree class function findRoot. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		
}

/*****************************************************************/
void Tree::printBranch(int node) {
	try {
		
		// you are not a leaf
		if (tree[node].getLChild() != -1) {
			out << "(";
			printBranch(tree[node].getLChild());
			out << ",";
			printBranch(tree[node].getRChild());
			out << ")";
		}else { //you are a leaf
			tree[node].printNode(out);  //prints out name and branch length
		}
		
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the Tree class Function printBranch. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the Tree class function printBranch. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		
}

/*****************************************************************/



