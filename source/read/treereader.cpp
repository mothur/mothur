//
//  treereader.cpp
//  Mothur
//
//  Created by Sarah Westcott on 4/11/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#include "treereader.h"
#include "readtree.h"
#include "groupmap.h"

/***********************************************************************/
TreeReader::TreeReader(string tf, string cf) : treefile(tf), countfile(cf)  { 
    try {
        m = MothurOut::getInstance();
        ct = new CountTable();
        ct->readTable(cf, true, false);
        Utils util;
        Treenames = util.parseTreeFile(treefile); //fills treenames
        
        //if no groupinfo in count file we need to add it
        if (!ct->hasGroupInfo()) {
            ct->addGroup("Group1");
            vector<string> namesOfSeqs = ct->getNamesOfSeqs();
            for (int i = 0; i < namesOfSeqs.size(); i++) { 
                ct->setAbund(namesOfSeqs[i], "Group1", ct->getNumSeqs(namesOfSeqs[i]));
            }
        }
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
TreeReader::TreeReader(string tf, string gf, string nf) : treefile(tf),  groupfile(gf), namefile(nf)  { 
    try {
        m = MothurOut::getInstance();
        Utils util;
        Treenames = util.parseTreeFile(treefile); //fills treenames
        countfile = "";
        ct = new CountTable();
        if (namefile != "") { ct->createTable(namefile, groupfile, true); }
        else {
            set<string> nameMap;
            map<string, string> groupMap;
            set<string> gps;
            for (int i = 0; i < Treenames.size(); i++) { nameMap.insert(Treenames[i]);  }
            if (groupfile == "") { gps.insert("Group1"); for (int i = 0; i < Treenames.size(); i++) { groupMap[Treenames[i]] = "Group1"; } }
            else {
                GroupMap g(groupfile); 
                g.readMap();
                vector<string> seqs = g.getNamesSeqs();
                for (int i = 0; i < seqs.size(); i++) {  
                    string group = g.getGroup(seqs[i]);
                    groupMap[seqs[i]] = group;
                    gps.insert(group);
                }
            }
            ct->createTable(nameMap, groupMap, gps);
        }

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
        
        int numUniquesInName = ct->getNumUniqueSeqs();
		//if (namefile != "") { numUniquesInName = readNamesFile(); }
		
		ReadTree* read = new ReadNewickTree(treefile, Treenames);
		int readOk = read->read(ct); 
		
		if (readOk != 0) { m->mothurOut("Read Terminated."); m->mothurOutEndLine();  delete read; m->setControl_pressed(true); return 0; }
		
		read->AssembleTrees();
		trees = read->getTrees();
		delete read;
        
		//make sure all files match
		//if you provide a namefile we will use the numNames in the namefile as long as the number of unique match the tree names size.
		int numNamesInTree;
		if (namefile != "")  {  
			if (numUniquesInName == Treenames.size()) {  numNamesInTree = ct->getNumSeqs();  }
			else {   numNamesInTree = Treenames.size();  }
		}else {  numNamesInTree = Treenames.size();  }
		
		
		//output any names that are in group file but not in tree
		if (numNamesInTree < ct->getNumSeqs()) {
            vector<string> namesSeqsCt = ct->getNamesOfSeqs();
			for (int i = 0; i < namesSeqsCt.size(); i++) {
				//is that name in the tree?
				int count = 0;
				for (int j = 0; j < Treenames.size(); j++) {
					if (namesSeqsCt[i] == Treenames[j]) { break; } //found it
					count++;
				}
				
				if (m->getControl_pressed()) { for (int i = 0; i < trees.size(); i++) { delete trees[i]; } return 0; }
				
				//then you did not find it so report it 
				if (count == Treenames.size()) {
                    m->mothurOut(namesSeqsCt[i] + " is in your name or group file and not in your tree. It will be disregarded."); m->mothurOutEndLine();
                    ct->remove(namesSeqsCt[i]);
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
/***********************************************************************/


