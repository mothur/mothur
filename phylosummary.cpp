/*
 *  rawTrainingDataMaker.cpp
 *  Mothur
 *
 *  Created by westcott on 4/21/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "phylosummary.h"

/**************************************************************************************************/

PhyloSummary::PhyloSummary(string refTfile, string groupFile){
	try {
		m = MothurOut::getInstance();
		maxLevel = 0;
		
		if (groupFile != "") {
			groupmap = new GroupMap(groupFile);
			groupmap->readMap();
		}else{
			groupmap = NULL;
		}
				
		//check for necessary files
		string taxFileNameTest = refTfile.substr(0,refTfile.find_last_of(".")+1) + "tree.sum";
		ifstream FileTest(taxFileNameTest.c_str());
		
		if (!FileTest) { 
			m->mothurOut("Error: can't find " + taxFileNameTest + "."); m->mothurOutEndLine(); exit(1);
		}else{
			readTreeStruct(FileTest);
		}
		
		tree[0].rank = "0";
		assignRank(0);

	}
	catch(exception& e) {
		m->errorOut(e, "PhyloSummary", "PhyloSummary");
		exit(1);
	}
}
/**************************************************************************************************/

void PhyloSummary::summarize(string userTfile){
	try {
		
		ifstream in;
		openInputFile(userTfile, in);
		
		//read in users taxonomy file and add sequences to tree
		string name, tax;
		while(!in.eof()){
			in >> name >> tax; gobble(in);
			
			addSeqToTree(name, tax);
			
			if (m->control_pressed) { break;  }
		}
		in.close();
	}
	catch(exception& e) {
		m->errorOut(e, "PhyloSummary", "summarize");
		exit(1);
	}
}

/**************************************************************************************************/

string PhyloSummary::getNextTaxon(string& heirarchy){
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
		m->errorOut(e, "PhyloSummary", "getNextTaxon");
		exit(1);
	}
}

/**************************************************************************************************/

int PhyloSummary::addSeqToTree(string seqName, string seqTaxonomy){
	try {
		numSeqs++;
		
		map<string, int>::iterator childPointer;
		
		int currentNode = 0;
		string taxon;
		
		int level = 0;
		
		while (seqTaxonomy != "") {
			
			if (m->control_pressed) { return 0; }
			
			//somehow the parent is getting one too many accnos
			//use print to reassign the taxa id
			taxon = getNextTaxon(seqTaxonomy);
			
			childPointer = tree[currentNode].children.find(taxon);
			
			if(childPointer != tree[currentNode].children.end()){	//if the node already exists, update count and move on
				if (groupmap != NULL) {
					//find out the sequences group
					string group = groupmap->getGroup(seqName);
					
					//do you have a count for this group?
					map<string, int>::iterator itGroup = tree[currentNode].groupCount.find(group);
					
					//if yes, increment it - there should not be a case where we can't find it since we load group in read
					if (itGroup != tree[currentNode].groupCount.end()) {
						tree[currentNode].groupCount[group]++;
					}
				}
				
				tree[currentNode].total++;

				currentNode = childPointer->second;
			}else{											//otherwise, create it
				m->mothurOut("Error: cannot find taxonomy in tree for " + seqName + "."); m->mothurOutEndLine();
				seqTaxonomy = "";
			}
			
			level++;
			
			if ((seqTaxonomy == "") && (level < maxLevel)) {  //if you think you are done and you are not.
				for (int k = level; k < maxLevel; k++) {  seqTaxonomy += "unclassified;";   }
			}
		}

	}
	catch(exception& e) {
		m->errorOut(e, "PhyloSummary", "addSeqToTree");
		exit(1);
	}
}
/**************************************************************************************************/

void PhyloSummary::assignRank(int index){
	try {
		map<string,int>::iterator it;
		int counter = 1;
		
		for(it=tree[index].children.begin();it!=tree[index].children.end();it++){
			tree[it->second].rank = tree[index].rank + '.' + toString(counter);
			counter++;
									
			assignRank(it->second);
		}
	}
	catch(exception& e) {
		m->errorOut(e, "PhyloSummary", "assignRank");
		exit(1);
	}
}
/**************************************************************************************************/

void PhyloSummary::print(ofstream& out){
	try {
		//print labels
		out << "taxlevel\t rank ID\t label\t daughterlevels\t total\t";
		if (groupmap != NULL) {
			for (int i = 0; i < groupmap->namesOfGroups.size(); i++) {
				out << groupmap->namesOfGroups[i] << '\t';
			}
		}
		
		out << endl;
		
		//print root
		out << tree[0].level << "\t" << tree[0].rank << "\t" << tree[0].name << "\t" << tree[0].children.size() << "\t" << tree[0].total << "\t";
		
		map<string, int>::iterator itGroup;
		if (groupmap != NULL) {
			for (itGroup = tree[0].groupCount.begin(); itGroup != tree[0].groupCount.end(); itGroup++) {
				out << itGroup->second << '\t';
			}
		}
		out << endl;
		
		//print rest
		print(0, out);
		
	}
	catch(exception& e) {
		m->errorOut(e, "PhyloSummary", "print");
		exit(1);
	}
}

/**************************************************************************************************/

void PhyloSummary::print(int i, ofstream& out){
	try {
		map<string,int>::iterator it;
		for(it=tree[i].children.begin();it!=tree[i].children.end();it++){
			
			if (tree[it->second].total != 0)  {
				out << tree[it->second].level << "\t" << tree[it->second].rank << "\t" << tree[it->second].name << "\t" << tree[it->second].children.size() << "\t" << tree[it->second].total << "\t";
				
				map<string, int>::iterator itGroup;
				if (groupmap != NULL) {
					for (itGroup = tree[it->second].groupCount.begin(); itGroup != tree[it->second].groupCount.end(); itGroup++) {
						out << itGroup->second << '\t';
					}
				}
				out << endl;
			}
			
			print(it->second, out);
		}
	}
	catch(exception& e) {
		m->errorOut(e, "PhyloSummary", "print");
		exit(1);
	}
}
/**************************************************************************************************/
void PhyloSummary::readTreeStruct(ifstream& in){
	try {
		int num;
		
		in >> num; gobble(in);
		
		tree.resize(num);
	
		//read the tree file
		for (int i = 0; i < tree.size(); i++) {
	
			in >> tree[i].level >> tree[i].name >> num; //num contains the number of children tree[i] has
			
			//set children
			string childName;
			int childIndex;
			for (int j = 0; j < num; j++) {
				in >> childName >> childIndex;
				tree[i].children[childName] = childIndex;
			}
			
			//initialize groupcounts
			if (groupmap != NULL) {
				for (int j = 0; j < groupmap->namesOfGroups.size(); j++) {
					tree[i].groupCount[groupmap->namesOfGroups[j]] = 0;
				}
			}
			
			tree[i].total = 0;
			
			gobble(in);
			
			if (tree[i].level > maxLevel) {  maxLevel = tree[i].level;  }
		}

	}
	catch(exception& e) {
		m->errorOut(e, "PhyloSummary", "print");
		exit(1);
	}
}

/**************************************************************************************************/


	
