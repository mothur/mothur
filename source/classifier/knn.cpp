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
Knn::Knn(string tfile, string tempFile, string method, int kmerSize, float gapOpen, float gapExtend, float match, float misMatch, int n, int tid) 
: Classify(), num(n), search(method) {
	try {
		threadID = tid;
        shortcuts = true;
		
		//create search database and names vector
		generateDatabaseAndNames(tfile, tempFile, method, kmerSize, gapOpen, gapExtend, match, misMatch);
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
		m->openOutputFile(outDistName, outDistance);
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
		 if (database != NULL) {  delete database; }
	}
	catch(exception& e) {
		m->errorOut(e, "Knn", "~Knn");
		exit(1);
	}
}
/**************************************************************************************************/
string Knn::getTaxonomy(Sequence* seq) {
	try {
		string tax;
		
		//use database to find closest seq
		vector<int> closest = database->findClosestSequences(seq, num);
	
		if (search == "distance") { ofstream outDistance; m->openOutputFileAppend(outDistName, outDistance); outDistance << seq->getName() << '\t' << database->getName(closest[0]) << '\t' << database->getSearchScore() << endl; outDistance.close();  }
	
		if (m->getControl_pressed()) { return tax; }

		vector<string> closestNames;
		for (int i = 0; i < closest.size(); i++) {
			//find that sequences taxonomy in map
			it = taxonomy.find(names[closest[i]]);
		
			//is this sequence in the taxonomy file
			if (it == taxonomy.end()) { //error not in file
				m->mothurOut("Error: sequence " + names[closest[i]] + " is not in the taxonomy file.  It will be eliminated as a match to sequence " + seq->getName() + "."); m->mothurOutEndLine();
			}else{   closestNames.push_back(it->first);	}
		}
		
		if (closestNames.size() == 0) {
			m->mothurOut("Error: All the matches for sequence " + seq->getName() + " have been eliminated. "); m->mothurOutEndLine();
			tax = "unknown;";
		}else{
			tax = findCommonTaxonomy(closestNames);
			if (tax == "") { m->mothurOut("There are no common levels for sequence " + seq->getName() + ". "); m->mothurOutEndLine(); tax = "unknown;"; }
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
		/*vector< vector<string> > taxons;  //taxon[0] = vector of taxonomy info for closest[0].
										//so if closest[0] taxonomy is Bacteria;Alphaproteobacteria;Rhizobiales;Azorhizobium_et_rel.;Methylobacterium_et_rel.;Bosea;
										//taxon[0][0] = Bacteria, taxon[0][1] = Alphaproteobacteria....
										
		taxons.resize(closest.size());
		int smallest = 100;
		
		for (int i = 0; i < closest.size(); i++) {
			if (m->getControl_pressed()) { return "control"; }
		
			string tax = taxonomy[closest[i]];  //we know its there since we checked in getTaxonomy
			cout << tax << endl;
		
			taxons[i] = parseTax(tax);
		
			//figure out who has the shortest taxonomy info. so you can start comparing there
			if (taxons[i].size() < smallest) {
				smallest = taxons[i].size();
			}
		}
	
		//start at the highest level all the closest seqs have
		string common = "";
		for (int i = (smallest-1); i >= 0; i--) {
			if (m->getControl_pressed()) { return "control"; }

			string thistax = taxons[0][i];
			int num = 0;
			for (int j = 1; j < taxons.size(); j++) {
				if (taxons[j][i] != thistax) { break; }
				num = j;
			}
		
			if (num == (taxons.size()-1)) { //they all match at this level
				for (int k = 0; k <= i; k++) {
					common += taxons[0][k] + ';';
				}
				break;
			}
		}*/
		
		string conTax;
		
		//create a tree containing sequences from this bin
		PhyloTree* p = new PhyloTree();
		
		for (int i = 0; i < closest.size(); i++) {
			p->addSeqToTree(closest[i], taxonomy[closest[i]]);
		}
		
		//build tree
		p->assignHeirarchyIDs(0);
		
		TaxNode currentNode = p->get(0);
		
		//at each level
		while (currentNode.children.size() != 0) { //you still have more to explore
			
			TaxNode bestChild;
			int bestChildSize = 0;
			
			//go through children
			for (map<string, int>::iterator itChild = currentNode.children.begin(); itChild != currentNode.children.end(); itChild++) {
				
				TaxNode temp = p->get(itChild->second);
				
				//select child with largest accessions - most seqs assigned to it
				if (temp.accessions.size() > bestChildSize) {
					bestChild = p->get(itChild->second);
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
		
		delete p;	
		
		return conTax;
	}
	catch(exception& e) {
		m->errorOut(e, "Knn", "findCommonTaxonomy");
		exit(1);
	}
}
/**************************************************************************************************/

