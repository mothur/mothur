//
//  treereader.cpp
//  Mothur
//
//  Created by Sarah Westcott on 4/11/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#include "treereader.h"
#include "readtree.h"

/***********************************************************************/

TreeReader::TreeReader(string tf) : treefile(tf)  { 
    try {
        m = MothurOut::getInstance();
        namefile = "";
        groupfile = "";
        readTrees();
    }
	catch(exception& e) {
		m->errorOut(e, "TreeReader", "TreeReader");
		exit(1);
	}
}
/***********************************************************************/

TreeReader::TreeReader(string tf, string gf) : treefile(tf),  groupfile(gf)  { 
    try {
        m = MothurOut::getInstance();
        namefile = "";
        readTrees();
    }
	catch(exception& e) {
		m->errorOut(e, "TreeReader", "TreeReader");
		exit(1);
	}
}
/***********************************************************************/
TreeReader::TreeReader(string tf, string gf, string nf) : treefile(tf),  groupfile(gf), namefile(nf)  { 
    try {
        m = MothurOut::getInstance();
        readTrees();
    }
	catch(exception& e) {
		m->errorOut(e, "TreeReader", "TreeReader");
		exit(1);
	}
}
/***********************************************************************/
bool TreeReader::readTrees()  { 
    try {
        
        tmap = new TreeMap();
        if (groupfile != "") {      tmap->readMap(groupfile);        }
		else{ //fake out by putting everyone in one group
			Tree* tree = new Tree(treefile); delete tree;  //extracts names from tree to make faked out groupmap
			for (int i = 0; i < m->Treenames.size(); i++) { tmap->addSeq(m->Treenames[i], "Group1"); }
		}
		
        int numUniquesInName = 0;
		if (namefile != "") { numUniquesInName = readNamesFile(); }
		
		ReadTree* read = new ReadNewickTree(treefile);
		int readOk = read->read(tmap); 
		
		if (readOk != 0) { m->mothurOut("Read Terminated."); m->mothurOutEndLine();  delete read; m->control_pressed=true; return 0; }
		
		read->AssembleTrees(names);
		trees = read->getTrees();
		delete read;
        
		//make sure all files match
		//if you provide a namefile we will use the numNames in the namefile as long as the number of unique match the tree names size.
		int numNamesInTree;
		if (namefile != "")  {  
			if (numUniquesInName == m->Treenames.size()) {  numNamesInTree = nameMap.size();  }
			else {   numNamesInTree = m->Treenames.size();  }
		}else {  numNamesInTree = m->Treenames.size();  }
		
		
		//output any names that are in group file but not in tree
		if (numNamesInTree < tmap->getNumSeqs()) {
			for (int i = 0; i < tmap->namesOfSeqs.size(); i++) {
				//is that name in the tree?
				int count = 0;
				for (int j = 0; j < m->Treenames.size(); j++) {
					if (tmap->namesOfSeqs[i] == m->Treenames[j]) { break; } //found it
					count++;
				}
				
				if (m->control_pressed) { for (int i = 0; i < trees.size(); i++) { delete trees[i]; } return 0; }
				
				//then you did not find it so report it 
				if (count == m->Treenames.size()) { 
					//if it is in your namefile then don't remove
					map<string, string>::iterator it = nameMap.find(tmap->namesOfSeqs[i]);
					
					if (it == nameMap.end()) {
						m->mothurOut(tmap->namesOfSeqs[i] + " is in your groupfile and not in your tree. It will be disregarded."); m->mothurOutEndLine();
						tmap->removeSeq(tmap->namesOfSeqs[i]);
						i--; //need this because removeSeq removes name from namesOfSeqs
					}
				}
			}
		}
        
        return true;
    }
	catch(exception& e) {
		m->errorOut(e, "TreeReader", "readTrees");
		exit(1);
	}
}
/*****************************************************************/
int TreeReader::readNamesFile() {
	try {
		nameMap.clear();
        names.clear();
		int numUniquesInName = 0;
		
		ifstream in;
		m->openInputFile(namefile, in);
		
		string first, second;
		map<string, string>::iterator itNames;
		
		while(!in.eof()) {
			in >> first >> second; m->gobble(in);
			
			numUniquesInName++;
			
			itNames = nameMap.find(first);
			if (itNames == nameMap.end()) {  
				names[first] = second; 
				
				//we need a list of names in your namefile to use above when removing extra seqs above so we don't remove them
				vector<string> dupNames;
				m->splitAtComma(second, dupNames);
				
				for (int i = 0; i < dupNames.size(); i++) {	
					nameMap[dupNames[i]] = first; 
					if ((groupfile == "") && (i != 0)) { tmap->addSeq(dupNames[i], "Group1"); } 
				}
			}else {  m->mothurOut(first + " has already been seen in namefile, disregarding names file."); m->mothurOutEndLine(); in.close(); nameMap.clear(); names.clear(); namefile = ""; return 1; }			
		}
		in.close();
		
		return numUniquesInName;
	}
	catch(exception& e) {
		m->errorOut(e, "TreeReader", "readNamesFile");
		exit(1);
	}
}
/***********************************************************************/


