/*
 *  rawTrainingDataMaker.cpp
 *  Mothur
 *
 *  Created by westcott on 4/21/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "rawtrainingdatamaker.h"

/**************************************************************************************************/

RawTrainingDataMaker::RawTrainingDataMaker(){
	try {
		m = MothurOut::getInstance();
		numNodes = 1;
		numSeqs = 0;
		tree.push_back(rawTaxNode("Root"));
		tree[0].rank = "Root";
		maxLevel = 0;
	}
	catch(exception& e) {
		m->errorOut(e, "RawTrainingDataMaker", "RawTrainingDataMaker");
		exit(1);
	}
}
/**************************************************************************************************/

RawTrainingDataMaker::RawTrainingDataMaker(string tfile){
	try {
		m = MothurOut::getInstance();
		numNodes = 1;
		numSeqs = 0;
		tree.push_back(rawTaxNode("Root"));
		tree[0].rank = "Root";
		maxLevel = 0;
		
		ifstream in;
		openInputFile(tfile, in);
		
		//read in users taxonomy file and add sequences to tree
		string name, tax;
		while(!in.eof()){
			in >> name >> tax; gobble(in);
			
			addSeqToTree(name, tax);
		}
		in.close();
	
		assignRank(0);
	}
	catch(exception& e) {
		m->errorOut(e, "RawTrainingDataMaker", "RawTrainingDataMaker");
		exit(1);
	}
}

/**************************************************************************************************/

string RawTrainingDataMaker::getNextTaxon(string& heirarchy){
	try {
		string currentLevel = "";
		if(heirarchy != ""){
			int pos = heirarchy.find_first_of(';');
			currentLevel=heirarchy.substr(0,pos);
			if (pos != (heirarchy.length()-1)) {  heirarchy=heirarchy.substr(pos+1);  }
			else { heirarchy = ""; }
		}
		return currentLevel;
	}
	catch(exception& e) {
		m->errorOut(e, "RawTrainingDataMaker", "getNextTaxon");
		exit(1);
	}
}

/**************************************************************************************************/

int RawTrainingDataMaker::addSeqToTree(string seqName, string seqTaxonomy){
	try {
		numSeqs++;
		
		map<string, int>::iterator childPointer;
		
		int currentNode = 0;
		string taxon;
		
		while(seqTaxonomy != ""){
			
			if (m->control_pressed) { return 0; }
			
			//somehow the parent is getting one too many accnos
			//use print to reassign the taxa id
			taxon = getNextTaxon(seqTaxonomy);
			
			childPointer = tree[currentNode].children.find(taxon);
			
			if(childPointer != tree[currentNode].children.end()){	//if the node already exists, move on
				currentNode = childPointer->second;
			}else{											//otherwise, create it
				tree.push_back(rawTaxNode(taxon));
				numNodes++;
				tree[currentNode].children[taxon] = numNodes-1;
				tree[numNodes-1].parent = currentNode;
				
				currentNode = tree[currentNode].children[taxon];
			}
		}

	}
	catch(exception& e) {
		m->errorOut(e, "RawTrainingDataMaker", "addSeqToTree");
		exit(1);
	}
}
/**************************************************************************************************/

void RawTrainingDataMaker::assignRank(int index){
	try {
		map<string,int>::iterator it;
				
		string ranks[9] = { "Root","Domain","Kingdom","Phylum","Class","Order","Family","Genus","Species" };
		
		for(it=tree[index].children.begin();it!=tree[index].children.end();it++){
			tree[it->second].level = tree[index].level + 1;
			
			if (tree[it->second].level > 8) { 
				tree[it->second].rank = ("unknown" + toString(tree[it->second].level));
			}else {
				tree[it->second].rank = ranks[tree[it->second].level];
			}
						
			//save maxLevel for binning the unclassified seqs
			if (tree[it->second].level > maxLevel) { maxLevel = tree[it->second].level; } 
			
			assignRank(it->second);
		}
	}
	catch(exception& e) {
		m->errorOut(e, "RawTrainingDataMaker", "assignRank");
		exit(1);
	}
}
/**************************************************************************************************/

void RawTrainingDataMaker::print(ofstream& out){
	try {
		//string temp = tree[0].name +" " + tree[0].rank;
		//sanityCheck[temp] = temp;
		
		out << "0" << "*" << tree[0].name << "*" << tree[0].parent << "*" << tree[0].level << "*" << tree[0].rank << endl;
		print(0, out);
		
	}
	catch(exception& e) {
		m->errorOut(e, "RawTrainingDataMaker", "print");
		exit(1);
	}
}

/**************************************************************************************************/

void RawTrainingDataMaker::print(int i, ofstream& out){
	try {
		map<string,int>::iterator it;
		for(it=tree[i].children.begin();it!=tree[i].children.end();it++){
			//string temp = tree[it->second].name + " " + tree[it->second].rank;
			
			//map<string, string>::iterator itSan;
			//itSan = sanityCheck.find(temp);
			
			//if (itSan == sanityCheck.end()) {
				out << it->second << "*" << tree[it->second].name << "*" << tree[it->second].parent << "*" << tree[it->second].level << "*" << tree[it->second].rank << endl;
				//sanityCheck[temp] = temp;
			//}
			print(it->second, out);
		}
	}
	catch(exception& e) {
		m->errorOut(e, "RawTrainingDataMaker", "print");
		exit(1);
	}
}
/**************************************************************************************************/


	
