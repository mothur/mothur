/*
 *  doTaxonomy.cpp
 *  
 *
 *  Created by Pat Schloss on 6/17/09.
 *  Copyright 2009 Patrick D. Schloss. All rights reserved.
 *
 */

#include "doTaxonomy.h"

/**************************************************************************************************/

PhyloTree::PhyloTree(){
	
	numNodes = 1;
	numSeqs = 0;
	tree.push_back(TaxNode("Root"));
	tree[0].heirarchyID = "0";
}

/**************************************************************************************************/

string PhyloTree::getNextTaxon(string& heirarchy){
	
	string currentLevel = "";
	if(heirarchy != ""){
		currentLevel=heirarchy.substr(0,heirarchy.find_first_of(';'));
		heirarchy=heirarchy.substr(heirarchy.find_first_of(';')+1);
	}
	return currentLevel;
	
}

/**************************************************************************************************/

void PhyloTree::addSeqToTree(string seqName, string seqTaxonomy){

	numSeqs++;

	map<string, int>::iterator childPointer;

	int currentNode = 0;
	int level = 1;
	
	tree[0].accessions.push_back(seqName);
	string taxon;// = getNextTaxon(seqTaxonomy);

	while(seqTaxonomy != ""){

		level++;

//somehow the parent is getting one too many accnos
//use print to reassign the taxa id
		taxon = getNextTaxon(seqTaxonomy);

		childPointer = tree[currentNode].children.find(taxon);

		if(childPointer != tree[currentNode].children.end()){	//if the node already exists, move on
			currentNode = childPointer->second;
			tree[currentNode].accessions.push_back(seqName);
		}
		else{											//otherwise, create it
			tree.push_back(TaxNode(taxon));
			numNodes++;
			tree[currentNode].children[taxon] = numNodes-1;

//			int numChildren = tree[currentNode].children.size();
//			string heirarchyID = tree[currentNode].heirarchyID;
//			tree[currentNode].accessions.push_back(seqName);
			
			currentNode = tree[currentNode].children[taxon];
			tree[currentNode].accessions.push_back(seqName);

//			tree[currentNode].level = level;
//			tree[currentNode].childNumber = numChildren;
//			tree[currentNode].heirarchyID = heirarchyID + '.' + toString(tree[currentNode].childNumber);
		}

	}
}

/**************************************************************************************************/

void PhyloTree::assignHeirarchyIDs(int index){
	
	map<string,int>::iterator it;
	int counter = 1;
	
	for(it=tree[index].children.begin();it!=tree[index].children.end();it++){
		tree[it->second].heirarchyID = tree[index].heirarchyID + '.' + toString(counter);
		counter++;
		tree[it->second].level = tree[index].level + 1;
		assignHeirarchyIDs(it->second);

	}
	
}

/**************************************************************************************************/

void PhyloTree::print(ofstream& out){
	
	out << tree[0].level << '\t'<< tree[0].heirarchyID << '\t' << tree[0].name << '\t' << tree[0].children.size() << '\t' << tree[0].accessions.size() << endl;
	print(0, out);

	
}

/**************************************************************************************************/

void PhyloTree::print(int i, ofstream& out){
	
	map<string,int>::iterator it;
	for(it=tree[i].children.begin();it!=tree[i].children.end();it++){
		out <<tree[it->second].level << '\t' << tree[it->second].heirarchyID << '\t' << tree[it->second].name << '\t' << tree[it->second].children.size() << '\t' << tree[it->second].accessions.size() << endl;
		print(it->second, out);
	}

}

/**************************************************************************************************/


	
