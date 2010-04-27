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
			
			iss >> numNodes; gobble(iss);
			
			tree.resize(numNodes);
			
			for (int i = 0; i < tree.size(); i++) {
				iss >> tree[i].name >> tree[i].level >> tree[i].parent; gobble(iss);
			}
			
			//read genus nodes
			int numGenus = 0;
			iss >> numGenus; gobble(iss);
			
			int gnode, gsize;
			totals.clear();
			for (int i = 0; i < numGenus; i++) {
				iss >> gnode >> gsize; gobble(iss);
				
				uniqueTaxonomies[gnode] = gnode;
				totals.push_back(gsize);
			}
			
			MPI_File_close(&inMPI);
			
		#else
			in >> numNodes; gobble(in);
			
			tree.resize(numNodes);
			
			for (int i = 0; i < tree.size(); i++) {
				in >> tree[i].name >> tree[i].level >> tree[i].parent; gobble(in);
			}
			
			//read genus nodes
			int numGenus = 0;
			in >> numGenus; gobble(in);
			
			int gnode, gsize;
			totals.clear();
			for (int i = 0; i < numGenus; i++) {
				in >> gnode >> gsize; gobble(in);
				
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
			int pid, num;
			vector<long> positions;
			
			MPI_Status status; 
			MPI_File inMPI;
			MPI_Comm_rank(MPI_COMM_WORLD, &pid); //find out who we are

			char inFileName[1024];
			strcpy(inFileName, tfile.c_str());

			MPI_File_open(MPI_COMM_WORLD, inFileName, MPI_MODE_RDONLY, MPI_INFO_NULL, &inMPI);  //comm, filename, mode, info, filepointer

			if (pid == 0) {
				positions = setFilePosEachLine(tfile, num);
				
				//send file positions to all processes
				MPI_Bcast(&num, 1, MPI_INT, 0, MPI_COMM_WORLD);  //send numSeqs
				MPI_Bcast(&positions[0], (num+1), MPI_LONG, 0, MPI_COMM_WORLD); //send file pos	
			}else{
				MPI_Bcast(&num, 1, MPI_INT, 0, MPI_COMM_WORLD); //get numSeqs
				positions.resize(num);
				MPI_Bcast(&positions[0], (num+1), MPI_LONG, 0, MPI_COMM_WORLD); //get file positions
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
		
		#else
			ifstream in;
			openInputFile(tfile, in);
			
			//read in users taxonomy file and add sequences to tree
			while(!in.eof()){
				in >> name >> tax; gobble(in);
				
				addSeqToTree(name, tax);
			}
			in.close();
		#endif
		
		assignHeirarchyIDs(0);
		
		//create file for summary if needed
		setUp(tfile);
	}
	catch(exception& e) {
		m->errorOut(e, "PhyloTree", "PhyloTree");
		exit(1);
	}
}

/**************************************************************************************************/

string PhyloTree::getNextTaxon(string& heirarchy){
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
		string taxon;// = getNextTaxon(seqTaxonomy);
		
		while(seqTaxonomy != ""){
			
			level++;
			
			if (m->control_pressed) { return 0; }
			
			//somehow the parent is getting one too many accnos
			//use print to reassign the taxa id
			taxon = getNextTaxon(seqTaxonomy);
			
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
				
				//			int numChildren = tree[currentNode].children.size();
				//			string heirarchyID = tree[currentNode].heirarchyID;
				//			tree[currentNode].accessions.push_back(seqName);
				
				currentNode = tree[currentNode].children[taxon];
				tree[currentNode].accessions.push_back(seqName);
				name2Taxonomy[seqName] = currentNode;
				//			tree[currentNode].level = level;
				//			tree[currentNode].childNumber = numChildren;
				//			tree[currentNode].heirarchyID = heirarchyID + '.' + toString(tree[currentNode].childNumber);
			}
		
			if (seqTaxonomy == "") {   uniqueTaxonomies[currentNode] = currentNode;	}
		}

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
			//create file needed for summary if it doesn't exist
			ifstream FileTest(taxFileNameTest.c_str());
			
			if (!FileTest) { 
				binUnclassified(taxFileNameTest); 
			}
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
		openOutputFile(file, out);
		
		map<string, int>::iterator itBin;
		map<string, int>::iterator childPointer;
		
		vector<TaxNode> copy = tree;
		int copyNodes = numNodes;
		
		//go through the seqs and if a sequence finest taxon is not the same level as the most finely defined taxon then classify it as unclassified where necessary
		for (itBin = name2Taxonomy.begin(); itBin != name2Taxonomy.end(); itBin++) {
			
			if (m->control_pressed) {  out.close(); break;  }
			
			int level = copy[itBin->second].level;
			int currentNode = itBin->second;
			
			//this sequence is unclassified at some levels
			while(level != maxLevel){
			
				level++;
			
				string taxon = "unclassified";	
				
				//does the parent have a child names 'unclassified'?
				childPointer = copy[currentNode].children.find(taxon);
				
				if(childPointer != copy[currentNode].children.end()){	//if the node already exists, move on
					currentNode = childPointer->second; //currentNode becomes 'unclassified'
					copy[currentNode].accessions.push_back(itBin->first);  //add this seq
				}
				else{											//otherwise, create it
					copy.push_back(TaxNode(taxon));
					copyNodes++;
					copy[currentNode].children[taxon] = copyNodes-1;
					copy[copyNodes-1].parent = currentNode;
					copy[copyNodes-1].level = copy[currentNode].level + 1;
									
					currentNode = copy[currentNode].children[taxon];
					copy[currentNode].accessions.push_back(itBin->first);
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
		out << copy.size() << endl;
		
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
			openOutputFile(treefilename, outTree);
			
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


	
