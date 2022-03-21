/*
 *  knn.cpp
 *  Mothur
 *
 *  Created by westcott on 11/4/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */

#include "knn.h"

/**************************************************************************************************/
Knn::Knn(string tfile, string tempFile, string method, int kmerSize, float gapOpen, float gapExtend, float match, float misMatch, int n, int tid, string version)
: Classify(), num(n), search(method) {
	try {
		threadID = tid;
        shortcuts = true;
		
		//create search database and names vector
		generateDatabaseAndNames(tfile, tempFile, method, kmerSize, gapOpen, gapExtend, match, misMatch, version);
	}
	catch(exception& e) {
		m->errorOut(e, "Knn", "Knn");
		exit(1);
	}
}
/**************************************************************************************************/
void Knn::setDistName(string s) {
	try {
		outDistName = s;
		ofstream outDistance;
        Utils util; util.openOutputFile(outDistName, outDistance);
		outDistance << "Name\tBestMatch\tDistance" << endl;
		outDistance.close();
	}
	catch(exception& e) {
		m->errorOut(e, "Knn", "setDistName");
		exit(1);
	}
}
/**************************************************************************************************/
Knn::~Knn() {
	try {
		 delete phyloTree; 
		 if (database != nullptr) {  delete database; }
	}
	catch(exception& e) {
		m->errorOut(e, "Knn", "~Knn");
		exit(1);
	}
}
/**************************************************************************************************/
string Knn::getTaxonomy(Sequence* seq, string& simpleTax, bool& flipped) {
	try {
		string tax;
		simpleTax = "";
		//use database to find closest seq
        vector<float> Scores;
		vector<int> closest = database->findClosestSequences(seq, num, Scores); 
	
         Utils util;
		if (search == "distance") { ofstream outDistance; util.openOutputFileAppend(outDistName, outDistance); outDistance << seq->getName() << '\t' << database->getName(closest[0]) << '\t' << Scores[0] << endl; outDistance.close();  }
	
		if (m->getControl_pressed()) { return tax; }

		vector<string> closestNames;
		for (int i = 0; i < closest.size(); i++) {
			//find that sequences taxonomy in map
			it = taxonomy.find(names[closest[i]]);
		
			//is this sequence in the taxonomy file
			if (it == taxonomy.end()) { //error not in file
				m->mothurOut("Error: sequence " + names[closest[i]] + " is not in the taxonomy file.  It will be eliminated as a match to sequence " + seq->getName() + ".\n"); 
			}else{   closestNames.push_back(it->first);	}
		}
		
		if (closestNames.size() == 0) {
			m->mothurOut("Error: All the matches for sequence " + seq->getName() + " have been eliminated. \n"); 
			tax = "unknown;";
		}else{
			tax = findCommonTaxonomy(closestNames);
			if (tax == "") { m->mothurOut("There are no common levels for sequence " + seq->getName() + ".\n"); tax = "unknown;"; }
		}
		
		simpleTax = tax;
		return tax;	
	}
	catch(exception& e) {
		m->errorOut(e, "Knn", "getTaxonomy");
		exit(1);
	}
}
/**************************************************************************************************/
string Knn::findCommonTaxonomy(vector<string> closest)  {
	try {
        string conTax;
		
		//create a tree containing sequences from this bin
		PhyloTree p;
		
		for (int i = 0; i < closest.size(); i++) { p.addSeqToTree(closest[i], taxonomy[closest[i]]); }
		
		//build tree
		p.assignHeirarchyIDs(0);
		
		TaxNode currentNode = p.get(0);
		
		//at each level
		while (currentNode.children.size() != 0) { //you still have more to explore
			
			TaxNode bestChild;
			int bestChildSize = 0;
			
			//go through children
			for (map<string, int>::iterator itChild = currentNode.children.begin(); itChild != currentNode.children.end(); itChild++) {
				
				TaxNode temp = p.get(itChild->second);
				
				//select child with largest accessions - most seqs assigned to it
				if (temp.accessions.size() > bestChildSize) {
					bestChild = p.get(itChild->second);
					bestChildSize = temp.accessions.size();
				}
				
			}
			
			if (bestChildSize == closest.size()) { //if yes, add it
				conTax += bestChild.name + ";";
			}else{ //if no, quit
				break;
			}
			
			//move down a level
			currentNode = bestChild;
		}
		
		return conTax;
	}
	catch(exception& e) {
		m->errorOut(e, "Knn", "findCommonTaxonomy");
		exit(1);
	}
}
/**************************************************************************************************/

