/*
 *  doTaxonomy.cpp
 *  
 *
 *  Created by Pat Schloss on 6/17/09.
 *  Copyright 2009 Patrick D. Schloss. All rights reserved.
 *
 */

#include "phylotree.h"

/**************************************************************************************************/

PhyloTree::PhyloTree(){
	try {
		m = MothurOut::getInstance();
		numNodes = 1;
		numSeqs = 0;
		tree.push_back(TaxNode("Root"));
		tree[0].heirarchyID = "0";
		maxLevel = 0;
		calcTotals = true;
		addSeqToTree("unknown", "unknown;");
	}
	catch(exception& e) {
		m->errorOut(e, "PhyloTree", "PhyloTree");
		exit(1);
	}
}
/**************************************************************************************************/

PhyloTree::PhyloTree(ifstream& in, string filename){
	try {
		m = MothurOut::getInstance();
		calcTotals = false;
		numNodes = 0;
		numSeqs = 0;
		
		#ifdef USE_MPI
			MPI_File inMPI;
			MPI_Offset size;
			MPI_Status status;

			char inFileName[1024];
			strcpy(inFileName, filename.c_str());

			MPI_File_open(MPI_COMM_WORLD, inFileName, MPI_MODE_RDONLY, MPI_INFO_NULL, &inMPI);  
			MPI_File_get_size(inMPI, &size);
			
			char* buffer = new char[size];
			MPI_File_read(inMPI, buffer, size, MPI_CHAR, &status);

			string tempBuf = buffer;
			if (tempBuf.length() > size) { tempBuf = tempBuf.substr(0, size);  }
			istringstream iss (tempBuf,istringstream::in);
			delete buffer;
			
			//read version
			m->getline(iss); m->gobble(iss);
			
			iss >> numNodes; m->gobble(iss);
			
			tree.resize(numNodes);
			
			for (int i = 0; i < tree.size(); i++) {
				iss >> tree[i].name >> tree[i].level >> tree[i].parent; m->gobble(iss);
			}
			
			//read genus nodes
			int numGenus = 0;
			iss >> numGenus; m->gobble(iss);
			
			int gnode, gsize;
			totals.clear();
			for (int i = 0; i < numGenus; i++) {
				iss >> gnode >> gsize; m->gobble(iss);
				
				uniqueTaxonomies[gnode] = gnode;
				totals.push_back(gsize);
			}
			
			MPI_File_close(&inMPI);
			
		#else
			//read version
			string line = m->getline(in); m->gobble(in);
			
			in >> numNodes; m->gobble(in);
			
			tree.resize(numNodes);
			
			for (int i = 0; i < tree.size(); i++) {
				in >> tree[i].name >> tree[i].level >> tree[i].parent; m->gobble(in);
			}
			
			//read genus nodes
			int numGenus = 0;
			in >> numGenus; m->gobble(in);
			
			int gnode, gsize;
			totals.clear();
			for (int i = 0; i < numGenus; i++) {
				in >> gnode >> gsize; m->gobble(in);
				
				uniqueTaxonomies[gnode] = gnode;
				totals.push_back(gsize);
			}
			
			in.close();
			
		#endif
		
	}
	catch(exception& e) {
		m->errorOut(e, "PhyloTree", "PhyloTree");
		exit(1);
	}
}
/**************************************************************************************************/

PhyloTree::PhyloTree(string tfile){
	try {
		m = MothurOut::getInstance();
		numNodes = 1;
		numSeqs = 0;
		tree.push_back(TaxNode("Root"));
		tree[0].heirarchyID = "0";
		maxLevel = 0;
		calcTotals = true;
		string name, tax;
		
		#ifdef USE_MPI
			int pid, num, processors;
			vector<unsigned long long> positions;
			
			MPI_Status status; 
			MPI_File inMPI;
			MPI_Comm_rank(MPI_COMM_WORLD, &pid); //find out who we are
			MPI_Comm_size(MPI_COMM_WORLD, &processors);

			char inFileName[1024];
			strcpy(inFileName, tfile.c_str());

			MPI_File_open(MPI_COMM_WORLD, inFileName, MPI_MODE_RDONLY, MPI_INFO_NULL, &inMPI);  //comm, filename, mode, info, filepointer

			if (pid == 0) {
				positions = m->setFilePosEachLine(tfile, num);
				
				//send file positions to all processes
				for(int i = 1; i < processors; i++) { 
					MPI_Send(&num, 1, MPI_INT, i, 2001, MPI_COMM_WORLD);
					MPI_Send(&positions[0], (num+1), MPI_LONG, i, 2001, MPI_COMM_WORLD);
				}
			}else{
				MPI_Recv(&num, 1, MPI_INT, 0, 2001, MPI_COMM_WORLD, &status);
				positions.resize(num+1);
				MPI_Recv(&positions[0], (num+1), MPI_LONG, 0, 2001, MPI_COMM_WORLD, &status);
			}
		
			//read file 
			for(int i=0;i<num;i++){
				//read next sequence
				int length = positions[i+1] - positions[i];
				char* buf4 = new char[length];

				MPI_File_read_at(inMPI, positions[i], buf4, length, MPI_CHAR, &status);

				string tempBuf = buf4;
				if (tempBuf.length() > length) { tempBuf = tempBuf.substr(0, length); }
				delete buf4;

				istringstream iss (tempBuf,istringstream::in);
				iss >> name >> tax;
				addSeqToTree(name, tax);
			}
			
			MPI_File_close(&inMPI);
			MPI_Barrier(MPI_COMM_WORLD); //make everyone wait - just in case
		
		#else
            map<string, string> temp;
            m->readTax(tfile, temp);
        
            for (map<string, string>::iterator itTemp = temp.begin(); itTemp != temp.end();) {
                addSeqToTree(itTemp->first, itTemp->second);
                temp.erase(itTemp++);
            }
		#endif
	
		assignHeirarchyIDs(0);
        
        
        string unknownTax = "unknown;";
        //added last taxon until you get desired level
		for (int i = 1; i < maxLevel; i++) {
			unknownTax += "unclassfied;";
		}
        
        addSeqToTree("unknown", unknownTax);
        
		//create file for summary if needed
		setUp(tfile);
	}
	catch(exception& e) {
		m->errorOut(e, "PhyloTree", "PhyloTree");
		exit(1);
	}
}

/**************************************************************************************************/

string PhyloTree::getNextTaxon(string& heirarchy, string seqname){
	try {
		string currentLevel = "";
		if(heirarchy != ""){
			int pos = heirarchy.find_first_of(';');
			
			if (pos == -1) { //you can't find another ;
				currentLevel = heirarchy;
				heirarchy = "";
				m->mothurOut(seqname + " is missing a ;, please check for other errors."); m->mothurOutEndLine();
			}else{
				currentLevel=heirarchy.substr(0,pos);
				if (pos != (heirarchy.length()-1)) {  heirarchy=heirarchy.substr(pos+1);  }
				else { heirarchy = ""; }
			}
			
		}
		return currentLevel;
	}
	catch(exception& e) {
		m->errorOut(e, "PhyloTree", "getNextTaxon");
		exit(1);
	}
}

/**************************************************************************************************/

int PhyloTree::addSeqToTree(string seqName, string seqTaxonomy){
	try {
		numSeqs++;
		
		map<string, int>::iterator childPointer;
		
		int currentNode = 0;
		int level = 1;
		
		tree[0].accessions.push_back(seqName);
		m->removeConfidences(seqTaxonomy);
		
		string taxon;// = getNextTaxon(seqTaxonomy);
	
		while(seqTaxonomy != ""){
			
			level++;
		
			if (m->control_pressed) { return 0; }
			
			//somehow the parent is getting one too many accnos
			//use print to reassign the taxa id
			taxon = getNextTaxon(seqTaxonomy, seqName);
			
			if (taxon == "") {  m->mothurOut(seqName + " has an error in the taxonomy.  This may be due to a ;;"); m->mothurOutEndLine(); if (currentNode != 0) {  uniqueTaxonomies[currentNode] = currentNode; } break;  }
			
			childPointer = tree[currentNode].children.find(taxon);
			
			if(childPointer != tree[currentNode].children.end()){	//if the node already exists, move on
				currentNode = childPointer->second;
				tree[currentNode].accessions.push_back(seqName);
				name2Taxonomy[seqName] = currentNode;
			}
			else{											//otherwise, create it
				tree.push_back(TaxNode(taxon));
				numNodes++;
				tree[currentNode].children[taxon] = numNodes-1;
				tree[numNodes-1].parent = currentNode;
				
				currentNode = tree[currentNode].children[taxon];
				tree[currentNode].accessions.push_back(seqName);
				name2Taxonomy[seqName] = currentNode;
			}
	
			if (seqTaxonomy == "") {   uniqueTaxonomies[currentNode] = currentNode;	}
		}
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "PhyloTree", "addSeqToTree");
		exit(1);
	}
}
/**************************************************************************************************/
vector<int> PhyloTree::getGenusNodes()	{
	try {
		genusIndex.clear();
		//generate genusIndexes
		map<int, int>::iterator it2;
		for (it2=uniqueTaxonomies.begin(); it2!=uniqueTaxonomies.end(); it2++) {  genusIndex.push_back(it2->first);	}
		
		return genusIndex;
	}
	catch(exception& e) {
		m->errorOut(e, "PhyloTree", "getGenusNodes");
		exit(1);
	}
}
/**************************************************************************************************/
vector<int> PhyloTree::getGenusTotals()	{
	try {
	
		if (calcTotals) {
			totals.clear();
			//reset counts because we are on a new word
			for (int j = 0; j < genusIndex.size(); j++) {
				totals.push_back(tree[genusIndex[j]].accessions.size());
			}
			return totals;
		}else{
			return totals;
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "PhyloTree", "getGenusNodes");
		exit(1);
	}
}
/**************************************************************************************************/

void PhyloTree::assignHeirarchyIDs(int index){
	try {
		map<string,int>::iterator it;
		int counter = 1;
		
		for(it=tree[index].children.begin();it!=tree[index].children.end();it++){
			tree[it->second].heirarchyID = tree[index].heirarchyID + '.' + toString(counter);
			counter++;
			tree[it->second].level = tree[index].level + 1;
						
			//save maxLevel for binning the unclassified seqs
			if (tree[it->second].level > maxLevel) { maxLevel = tree[it->second].level; } 
			
			assignHeirarchyIDs(it->second);
		}
	}
	catch(exception& e) {
		m->errorOut(e, "PhyloTree", "assignHeirarchyIDs");
		exit(1);
	}
}
/**************************************************************************************************/
void PhyloTree::setUp(string tfile){
	try{
		string taxFileNameTest = tfile.substr(0,tfile.find_last_of(".")+1) + "tree.sum";
		
		#ifdef USE_MPI
			int pid;
			MPI_Comm_rank(MPI_COMM_WORLD, &pid); //find out who we are

			if (pid == 0) {  binUnclassified(taxFileNameTest);  }
		
		#else
			binUnclassified(taxFileNameTest); 
		#endif
	}
	catch(exception& e) {
		m->errorOut(e, "PhyloTree", "setUp");
		exit(1);
	}
}
/**************************************************************************************************/
void PhyloTree::binUnclassified(string file){
	try {
	
		ofstream out;
		m->openOutputFile(file, out);
		
		map<string, int>::iterator itBin;
		map<string, int>::iterator childPointer;
		
		vector<TaxNode> copy = tree;
		
		//fill out tree
		fillOutTree(0, copy);
	
		//get leaf nodes that may need extension
		for (int i = 0; i < copy.size(); i++) {  

			if (copy[i].children.size() == 0) {
				leafNodes[i] = i;
			}
		}
		
		int copyNodes = copy.size();
	
		//go through the seqs and if a sequence finest taxon is not the same level as the most finely defined taxon then classify it as unclassified where necessary
		map<int, int>::iterator itLeaf;
		for (itLeaf = leafNodes.begin(); itLeaf != leafNodes.end(); itLeaf++) {
			
			if (m->control_pressed) {  out.close(); break;  }
			
			int level = copy[itLeaf->second].level;
			int currentNode = itLeaf->second;
			
			//this sequence is unclassified at some levels
			while(level < maxLevel){
		
				level++;
			
				string taxon = "unclassified";	
				
				//does the parent have a child names 'unclassified'?
				childPointer = copy[currentNode].children.find(taxon);
				
				if(childPointer != copy[currentNode].children.end()){	//if the node already exists, move on
					currentNode = childPointer->second; //currentNode becomes 'unclassified'
				}
				else{											//otherwise, create it
					copy.push_back(TaxNode(taxon));
					copyNodes++;
					copy[currentNode].children[taxon] = copyNodes-1;
					copy[copyNodes-1].parent = currentNode;
					copy[copyNodes-1].level = copy[currentNode].level + 1;
									
					currentNode = copy[currentNode].children[taxon];
				}
			}
		}
		
		if (!m->control_pressed) {
			//print copy tree
			print(out, copy);
		}
				
	}
	catch(exception& e) {
		m->errorOut(e, "PhyloTree", "binUnclassified");
		exit(1);
	}
}
/**************************************************************************************************/
void PhyloTree::fillOutTree(int index, vector<TaxNode>& copy) {
	try {
	
		map<string,int>::iterator it;
		
		it = copy[index].children.find("unclassified");
		if (it == copy[index].children.end()) { //no unclassified at this level
			string taxon = "unclassified";
			copy.push_back(TaxNode(taxon));
			copy[index].children[taxon] = copy.size()-1;
			copy[copy.size()-1].parent = index;
			copy[copy.size()-1].level = copy[index].level + 1;
		}
		
		if (tree[index].level < maxLevel) {
			for(it=tree[index].children.begin();it!=tree[index].children.end();it++){ //check your children
				fillOutTree(it->second, copy);
			}
		}

	}
	catch(exception& e) {
		m->errorOut(e, "PhyloTree", "fillOutTree");
		exit(1);
	}
}
/**************************************************************************************************/
string PhyloTree::getFullTaxonomy(string seqName) {
	try {
		string tax = "";
		
		int currentNode = name2Taxonomy[seqName];
		
		while (tree[currentNode].parent != -1) {
			tax = tree[currentNode].name + ";" + tax;
			currentNode = tree[currentNode].parent;
		}
		
		return tax;
	}
	catch(exception& e) {
		m->errorOut(e, "PhyloTree", "getFullTaxonomy");
		exit(1);
	}
}
/**************************************************************************************************/

void PhyloTree::print(ofstream& out, vector<TaxNode>& copy){
	try {
		
		//output mothur version
		out << "#" << m->getVersion() << endl;
		
		out << copy.size() << endl;
		
		out << maxLevel << endl;
				
		for (int i = 0; i < copy.size(); i++) {
				
			out << copy[i].level << '\t'<< copy[i].name << '\t' << copy[i].children.size() << '\t';
			
			map<string,int>::iterator it;
			for(it=copy[i].children.begin();it!=copy[i].children.end();it++){
				out << it->first << '\t' << it->second << '\t';
			}
			out << endl;
		}
		
		out.close();
	}
	catch(exception& e) {
		m->errorOut(e, "PhyloTree", "print");
		exit(1);
	}
}
/**************************************************************************************************/
void PhyloTree::printTreeNodes(string treefilename) {
	try {
	
		#ifdef USE_MPI
			int pid;
			MPI_Comm_rank(MPI_COMM_WORLD, &pid); //find out who we are

			if (pid == 0) {  
		
		#endif

			ofstream outTree;
			m->openOutputFile(treefilename, outTree);
			
			//output mothur version
			outTree << "#" << m->getVersion() << endl;
			
			//print treenodes
			outTree << tree.size() << endl;
			for (int i = 0; i < tree.size(); i++) {
				outTree << tree[i].name << '\t' << tree[i].level << '\t' << tree[i].parent << endl;
			}
			
			//print genus nodes
			outTree << endl << uniqueTaxonomies.size() << endl;
			map<int, int>::iterator it2;
			for (it2=uniqueTaxonomies.begin(); it2!=uniqueTaxonomies.end(); it2++) {  outTree << it2->first << '\t' << tree[it2->first].accessions.size() << endl;	}
			outTree << endl;
			
			outTree.close();
		
		#ifdef USE_MPI
			}
		#endif

		
	}
	catch(exception& e) {
		m->errorOut(e, "PhyloTree", "printTreeNodes");
		exit(1);
	}
}
/**************************************************************************************************/
TaxNode PhyloTree::get(int i ){
	try {
		if (i < tree.size()) {  return tree[i];	 }
		else {  cout << i << '\t' << tree.size() << endl ; m->mothurOut("Mismatch with taxonomy and template files. Cannot continue."); m->mothurOutEndLine(); exit(1); }
	}
	catch(exception& e) {
		m->errorOut(e, "PhyloTree", "get");
		exit(1);
	}
}
/**************************************************************************************************/
TaxNode PhyloTree::get(string seqName){
	try {
		map<string, int>::iterator itFind = name2Taxonomy.find(seqName);
	
		if (itFind != name2Taxonomy.end()) {  return tree[name2Taxonomy[seqName]];  }
		else { m->mothurOut("Cannot find " + seqName + ". Mismatch with taxonomy and template files. Cannot continue."); m->mothurOutEndLine(); exit(1);}
	}
	catch(exception& e) {
		m->errorOut(e, "PhyloTree", "get");
		exit(1);
	}
}
/**************************************************************************************************/
string PhyloTree::getName(int i ){
	try {
		if (i < tree.size()) {  return tree[i].name;	 }
		else { m->mothurOut("Mismatch with taxonomy and template files. Cannot continue."); m->mothurOutEndLine(); exit(1); }
	}
	catch(exception& e) {
		m->errorOut(e, "PhyloTree", "get");
		exit(1);
	}
}
/**************************************************************************************************/
int PhyloTree::getIndex(string seqName){
	try {
		map<string, int>::iterator itFind = name2Taxonomy.find(seqName);
	
		if (itFind != name2Taxonomy.end()) {  return name2Taxonomy[seqName];  }
		else { m->mothurOut("Cannot find " + seqName + ". Mismatch with taxonomy and template files. Cannot continue."); m->mothurOutEndLine(); exit(1);}
	}
	catch(exception& e) {
		m->errorOut(e, "PhyloTree", "get");
		exit(1);
	}
}
/**************************************************************************************************/
bool PhyloTree::ErrorCheck(vector<string> templateFileNames){
	try {
	
		bool okay = true;
		templateFileNames.push_back("unknown");
		
		map<string, int>::iterator itFind;
		map<string, int> taxonomyFileNames = name2Taxonomy;
		
        if (m->debug) { m->mothurOut("[DEBUG]: in error check. Numseqs in template = " + toString(templateFileNames.size()) + ". Numseqs in taxonomy = " + toString(taxonomyFileNames.size()) + ".\n"); }
        
		for (int i = 0; i < templateFileNames.size(); i++) {
			itFind = taxonomyFileNames.find(templateFileNames[i]);
			
			if (itFind != taxonomyFileNames.end()) { //found it so erase it
				taxonomyFileNames.erase(itFind);
			}else {
				m->mothurOut("'" +templateFileNames[i] + "' is in your template file and is not in your taxonomy file. Please correct."); m->mothurOutEndLine();
				okay = false;
			}
			
			//templateFileNames.erase(templateFileNames.begin()+i);
			//i--;
		}
		templateFileNames.clear();
		
		if (taxonomyFileNames.size() > 0) { //there are names in tax file that are not in template
			okay = false;
			
			for (itFind = taxonomyFileNames.begin(); itFind != taxonomyFileNames.end(); itFind++) {
				m->mothurOut(itFind->first + " is in your taxonomy file and is not in your template file. Please correct."); m->mothurOutEndLine();
			}
		}
		
		return okay;
	}
	catch(exception& e) {
		m->errorOut(e, "PhyloTree", "ErrorCheck");
		exit(1);
	}
}
/**************************************************************************************************/
	


	
