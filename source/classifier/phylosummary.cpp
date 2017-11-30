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

PhyloSummary::PhyloSummary(string refTfile, CountTable* c, bool r, int p){
	try {
		m = MothurOut::getInstance();
		maxLevel = 0;
		ignore = false;
        numSeqs = 0;
        relabund = r;
        printlevel = p;
		
		ct = c;
        groupmap = NULL;
        
		//check for necessary files
		string taxFileNameTest = util.getFullPathName((refTfile.substr(0,refTfile.find_last_of(".")+1) + "tree.sum"));
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

PhyloSummary::PhyloSummary(CountTable* c, bool r, int p){
	try {
		m = MothurOut::getInstance();
		maxLevel = 0;
		ignore = true;
        numSeqs = 0;
        relabund = r;
        printlevel = p;
		
		ct = c;
        groupmap = NULL;
		
		tree.push_back(rawTaxNode("Root"));
		tree[0].rank = "0";
	}
	catch(exception& e) {
		m->errorOut(e, "PhyloSummary", "PhyloSummary");
		exit(1);
	}
}
/**************************************************************************************************/
PhyloSummary::PhyloSummary(string refTfile, GroupMap* g, bool r, int p){
	try {
		m = MothurOut::getInstance();
		maxLevel = 0;
		ignore = false;
        numSeqs = 0;
        relabund = r;
        printlevel = p;
		
		groupmap = g;
        ct = NULL;
				
		//check for necessary files
		string taxFileNameTest = util.getFullPathName((refTfile.substr(0,refTfile.find_last_of(".")+1) + "tree.sum"));
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

PhyloSummary::PhyloSummary(GroupMap* g, bool r, int p){
	try {
		m = MothurOut::getInstance();
		maxLevel = 0;
		ignore = true;
        numSeqs = 0;
        relabund = r;
        printlevel = p;
		
		groupmap = g;
        ct = NULL;
		
		tree.push_back(rawTaxNode("Root"));
		tree[0].rank = "0";
	}
	catch(exception& e) {
		m->errorOut(e, "PhyloSummary", "PhyloSummary");
		exit(1);
	}
}
/**************************************************************************************************/

int PhyloSummary::summarize(string userTfile){
	try {
		map<string, string> temp;
        util.readTax(userTfile, temp, true);
        
        for (map<string, string>::iterator itTemp = temp.begin(); itTemp != temp.end();) {
            addSeqToTree(itTemp->first, itTemp->second);
            temp.erase(itTemp++);
        }
        
        return numSeqs;
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
		
		//are there confidence scores, if so remove them
		if (seqTaxonomy.find_first_of('(') != -1) {  util.removeConfidences(seqTaxonomy);	}
		
		while (seqTaxonomy != "") {
			
            level++;
            
			if (m->getControl_pressed()) { return 0; }
			
			//somehow the parent is getting one too many accnos
			//use print to reassign the taxa id
			taxon = getNextTaxon(seqTaxonomy);
			
			childPointer = tree[currentNode].children.find(taxon);
			
			if(childPointer != tree[currentNode].children.end()){	//if the node already exists, update count and move on
				int thisCount = 1;
                
                if (groupmap != NULL) {
					//find out the sequences group
					string group = groupmap->getGroup(seqName);
					
					if (group == "not found") {  m->mothurOut("[WARNING]: " + seqName + " is not in your groupfile, and will be included in the overall total, but not any group total."); m->mothurOutEndLine();  }
					
					//do you have a count for this group?
					map<string, int>::iterator itGroup = tree[childPointer->second].groupCount.find(group);
					
					//if yes, increment it - there should not be a case where we can't find it since we load group in read
					if (itGroup != tree[childPointer->second].groupCount.end()) {
						tree[childPointer->second].groupCount[group]++;
					}
				}else if (ct != NULL) {
                    if (ct->hasGroupInfo()) {
                        vector<int> groupCounts = ct->getGroupCounts(seqName);
                        vector<string> groups = ct->getNamesOfGroups();
                        for (int i = 0; i < groups.size(); i++) {
                            
                            if (groupCounts[i] != 0) {
                                //do you have a count for this group?
                                map<string, int>::iterator itGroup = tree[childPointer->second].groupCount.find(groups[i]);
                                
                                //if yes, increment it - there should not be a case where we can't find it since we load group in read
                                if (itGroup != tree[childPointer->second].groupCount.end()) {
                                    tree[childPointer->second].groupCount[groups[i]] += groupCounts[i];
                                }
                            }
                        }
                    }
                    thisCount = ct->getNumSeqs(seqName);
                }
				
				tree[childPointer->second].total += thisCount;

				currentNode = childPointer->second;
			}else{	
				if (ignore) {
						
					tree.push_back(rawTaxNode(taxon));
					int index = tree.size() - 1;
				
					tree[index].parent = currentNode;
					tree[index].level = level;
					tree[currentNode].children[taxon] = index;
                    int thisCount = 1;
					
					//initialize groupcounts
					if (groupmap != NULL) {
						vector<string> mGroups = groupmap->getNamesOfGroups();
						for (int j = 0; j < mGroups.size(); j++) {
							tree[index].groupCount[mGroups[j]] = 0;
						}
						
						//find out the sequences group
						string group = groupmap->getGroup(seqName);
						
						if (group == "not found") {  m->mothurOut("[WARNING]: " + seqName + " is not in your groupfile, and will be included in the overall total, but not any group total."); m->mothurOutEndLine();  }
						
						//do you have a count for this group?
						map<string, int>::iterator itGroup = tree[index].groupCount.find(group);
						
						//if yes, increment it - there should not be a case where we can't find it since we load group in read
						if (itGroup != tree[index].groupCount.end()) {
							tree[index].groupCount[group]++;
						}
					}else if (ct != NULL) {
                        if (ct->hasGroupInfo()) {
                            vector<string> mGroups = ct->getNamesOfGroups();
                            for (int j = 0; j < mGroups.size(); j++) {
                                tree[index].groupCount[mGroups[j]] = 0;
                            }
                            vector<int> groupCounts = ct->getGroupCounts(seqName);
                            vector<string> groups = ct->getNamesOfGroups();
                        
                            for (int i = 0; i < groups.size(); i++) {
                                if (groupCounts[i] != 0) {
                                   
                                    //do you have a count for this group?
                                    map<string, int>::iterator itGroup = tree[index].groupCount.find(groups[i]);
                                     
                                    //if yes, increment it - there should not be a case where we can't find it since we load group in read
                                    if (itGroup != tree[index].groupCount.end()) {
                                        tree[index].groupCount[groups[i]]+=groupCounts[i];
                                    }
                                }
                            }
                        }
                        thisCount = ct->getNumSeqs(seqName);
                    }
					
                    tree[index].total = thisCount;
					currentNode = index;
					
				}else{ //otherwise, error
					m->mothurOut("Warning: cannot find taxon " + taxon + " in reference taxonomy tree at level " + toString(tree[currentNode].level) + " for " + seqName + ". This may cause totals of daughter levels not to add up in summary file."); m->mothurOutEndLine();
					break;
				}
			}
        }
        
        if (level > maxLevel) { maxLevel = level; }
        
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "PhyloSummary", "addSeqToTree");
		exit(1);
	}
}
/**************************************************************************************************/

int PhyloSummary::addSeqToTree(string seqTaxonomy, map<string, bool> containsGroup){
	try {
		numSeqs++;
		
		map<string, int>::iterator childPointer;
		
		int currentNode = 0;
		string taxon;
		
		int level = 0;
		
		//are there confidence scores, if so remove them
		if (seqTaxonomy.find_first_of('(') != -1) {  util.removeConfidences(seqTaxonomy);	}
		
		while (seqTaxonomy != "") {
			
            level++;
            
			if (m->getControl_pressed()) { return 0; }
			
			//somehow the parent is getting one too many accnos
			//use print to reassign the taxa id
			taxon = getNextTaxon(seqTaxonomy);
			
			childPointer = tree[currentNode].children.find(taxon);
			
			if(childPointer != tree[currentNode].children.end()){	//if the node already exists, update count and move on
                for (map<string, bool>::iterator itGroup = containsGroup.begin(); itGroup != containsGroup.end(); itGroup++) {
                    if (itGroup->second ) {
                        tree[childPointer->second].groupCount[itGroup->first]++;
                    }
                }
					
				tree[childPointer->second].total++;
				
				currentNode = childPointer->second;
			}else{	
				if (ignore) {
					
					tree.push_back(rawTaxNode(taxon));
					int index = tree.size() - 1;
					
					tree[index].parent = currentNode;
					tree[index].level = level;
					tree[index].total = 1;
					tree[currentNode].children[taxon] = index;
						
                    for (map<string, bool>::iterator itGroup = containsGroup.begin(); itGroup != containsGroup.end(); itGroup++) {
                        if (itGroup->second ) {
                            tree[index].groupCount[itGroup->first]++;
                        }
                    }
					
					currentNode = index;
					
				}else{ //otherwise, error
					m->mothurOut("Warning: cannot find taxon " + taxon + " in reference taxonomy tree at level " + toString(tree[currentNode].level) + ". This may cause totals of daughter levels not to add up in summary file."); m->mothurOutEndLine();
					break;
				}
			}
		}
        
        if (level > maxLevel) { maxLevel = level; }
        
		return 0;
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

string PhyloSummary::findTaxon(string rank){
    try {
        vector<string> pieces; vector<int> indexes;
        util.splitAtChar(rank, pieces, '.');
        for (int i = 0; i < pieces.size(); i++) {
            int temp;
            util.mothurConvert(pieces[i], temp);
            indexes.push_back(temp);
        }
        string taxon = "";
        getTaxons(indexes, 1, 0, taxon);
        
        return taxon;
    }
    catch(exception& e) {
        m->errorOut(e, "PhyloSummary", "findTaxon");
        exit(1);
    }
}
/**************************************************************************************************/

string PhyloSummary::getTaxons(vector<int> indexes, int index, int i, string& taxon){
    try {
        int counter = 1;
        
        for(map<string,int>::iterator it=tree[i].children.begin();it!=tree[i].children.end();it++){
            if (counter == indexes[index]) {
                taxon += tree[it->second].name + ";";
                getTaxons(indexes, index+1, it->second, taxon);
            }
            counter++;
        }
        
        return taxon;
    }
    catch(exception& e) {
        m->errorOut(e, "PhyloSummary", "getNextTaxon");
        exit(1);
    }
}

/**************************************************************************************************/

void PhyloSummary::print(ofstream& out, string output){
	try {
		
		if (ignore)     {  assignRank(0); }
        vector<string> mGroups;
        
        //print labels
        if (output == "detail") {   out << "taxlevel\trankID\ttaxon\tdaughterlevels\ttotal";  }
        else                    {   out << "taxonomy\ttotal";  }
        
        if (printlevel == -1) { printlevel = maxLevel; }
        else if (printlevel > maxLevel) { m->mothurOut("[WARNING]: Your printlevel is greater than your maxlevel, adjusting your printlevel to " + toString(maxLevel) + "\n"); printlevel = maxLevel; }
        
		if (groupmap != NULL) {
			//so the labels match the counts below, since the map sorts them automatically...
			//sort(groupmap->namesOfGroups.begin(), groupmap->namesOfGroups.end());
            mGroups = groupmap->getNamesOfGroups();
			for (int i = 0; i < mGroups.size(); i++) {
				out << '\t' << mGroups[i];
			}
		}else if (ct != NULL) {
            if (ct->hasGroupInfo()) {
                mGroups = ct->getNamesOfGroups();
                for (int i = 0; i < mGroups.size(); i++) {
                    out << '\t' << mGroups[i];
                }
            }
        }
		out << endl;
		
		int totalChildrenInTree = 0;
		map<string, int>::iterator itGroup;
		map<string,int>::iterator it;
		for(it=tree[0].children.begin();it!=tree[0].children.end();it++){   
			if (tree[it->second].total != 0)  {   
				totalChildrenInTree++; 
				tree[0].total += tree[it->second].total;
				
				if (groupmap != NULL) {
					for (int i = 0; i < mGroups.size(); i++) { tree[0].groupCount[mGroups[i]] += tree[it->second].groupCount[mGroups[i]]; } 
				}else if ( ct != NULL) {
                    if (ct->hasGroupInfo()) { for (int i = 0; i < mGroups.size(); i++) { tree[0].groupCount[mGroups[i]] += tree[it->second].groupCount[mGroups[i]]; } }
                }
			}
		}
		
        
        
            //print root
            if (relabund) {
                out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);
                
                if (output == "detail") {
                    out << tree[0].level << "\t" << tree[0].rank << "\t" << tree[0].name << "\t" << totalChildrenInTree << "\t" << (tree[0].total/(double) tree[0].total);
                }else{
                    out << tree[0].name << "\t" << (tree[0].total/(double) tree[0].total);
                }
                
                
                if (groupmap != NULL) {
                    for (int i = 0; i < mGroups.size(); i++) {
                        double thisNum = tree[0].groupCount[mGroups[i]];
                        thisNum /= (double) groupmap->getNumSeqs(mGroups[i]);
                        out  << '\t' << thisNum;
                    }
                }else if ( ct != NULL) {
                    if (ct->hasGroupInfo()) {
                        for (int i = 0; i < mGroups.size(); i++) {
                            double thisNum = tree[0].groupCount[mGroups[i]];
                            thisNum /= (double) ct->getGroupCount(mGroups[i]);
                            out  << '\t' << thisNum;
                        }
                    }
                }
                out << endl;
               
                
            }else {
                if (output == "detail") {
                    out << tree[0].level << "\t" << tree[0].rank << "\t" << tree[0].name << "\t" << totalChildrenInTree << "\t" << tree[0].total;
                }else{
                    out << tree[0].name << '\t' << tree[0].total;
                }
                
                if (groupmap != NULL) {
                    for (int i = 0; i < mGroups.size(); i++) {  out  << '\t'<< tree[0].groupCount[mGroups[i]]; }
                }else if ( ct != NULL) {
                    if (ct->hasGroupInfo()) { for (int i = 0; i < mGroups.size(); i++) {  out  << '\t' << tree[0].groupCount[mGroups[i]]; } }
                }
                out << endl;
                
            }
        
        
		//print rest
		print(0, out, output);
		
	}
	catch(exception& e) {
		m->errorOut(e, "PhyloSummary", "print");
		exit(1);
	}
}
/**************************************************************************************************/

void PhyloSummary::print(ofstream& out, bool relabund){
	try {
		
		if (ignore) { assignRank(0); }
	
		int totalChildrenInTree = 0;
		map<string, int>::iterator itGroup;
		
		map<string,int>::iterator it;
		for(it=tree[0].children.begin();it!=tree[0].children.end();it++){
			if (tree[it->second].total != 0)  {
				totalChildrenInTree++;
				tree[0].total += tree[it->second].total;
				
				if (groupmap != NULL) {
                    vector<string> mGroups = groupmap->getNamesOfGroups();
					for (int i = 0; i < mGroups.size(); i++) { tree[0].groupCount[mGroups[i]] += tree[it->second].groupCount[mGroups[i]]; }
				}else if ( ct != NULL) {
                    vector<string> mGroups = ct->getNamesOfGroups();
                    if (ct->hasGroupInfo()) { for (int i = 0; i < mGroups.size(); i++) { tree[0].groupCount[mGroups[i]] += tree[it->second].groupCount[mGroups[i]]; } }
                }
			}
		}
        
        //print root
        out << tree[0].name << "\t" << "1.0000"; //root relative abundance is 1, everyone classifies to root
        
        if (groupmap != NULL) {
            vector<string> mGroups = groupmap->getNamesOfGroups();
            for (int i = 0; i < mGroups.size(); i++) {  out  << '\t' << "1.0000"; }
        }else if ( ct != NULL) {
            vector<string> mGroups = ct->getNamesOfGroups();
            if (ct->hasGroupInfo()) { for (int i = 0; i < mGroups.size(); i++) {  out  << '\t' << "1.0000"; } }
        }
        
        out << endl;
        
       
        
		//print rest
		print(0, out, relabund);
		
	}
	catch(exception& e) {
		m->errorOut(e, "PhyloSummary", "print");
		exit(1);
	}
}
/**************************************************************************************************/

void PhyloSummary::print(int i, ofstream& out, string output){
	try {
		map<string,int>::iterator it;
		for(it=tree[i].children.begin();it!=tree[i].children.end();it++){
			
			if (tree[it->second].total != 0)  {
                
				int totalChildrenInTree = 0;
                
				map<string,int>::iterator it2;
				for(it2=tree[it->second].children.begin();it2!=tree[it->second].children.end();it2++){
					if (tree[it2->second].total != 0)  {   totalChildrenInTree++; }
				}
                
                if ((output == "detail") && (printlevel >= tree[it->second].level)) {
                    if (relabund) {
                        out << tree[it->second].level << "\t" << tree[it->second].rank << "\t" << tree[it->second].name << "\t" << totalChildrenInTree << "\t" << (tree[it->second].total/(double) tree[0].total);
                    }else {
                        out << tree[it->second].level << "\t" << tree[it->second].rank << "\t" << tree[it->second].name << "\t" << totalChildrenInTree << "\t" << tree[it->second].total;
                    }
                }else {
                    if (printlevel == tree[it->second].level) { //leaf node - we want to print it. Use rank to find full taxonomy
                        if (relabund) {
                            out << findTaxon(tree[it->second].rank) << '\t' << tree[it->second].total/(double) tree[0].total;
                        }else {
                            out << findTaxon(tree[it->second].rank) << '\t' << tree[it->second].total;
                        }
                    }
                }
                
                if (relabund) {
                    map<string, int>::iterator itGroup;
                    if (groupmap != NULL) {
                        vector<string> mGroups = groupmap->getNamesOfGroups();
                        for (int i = 0; i < mGroups.size(); i++) {
                            if (((output == "detail") && (printlevel >= tree[it->second].level)) || (printlevel == tree[it->second].level)) {
                                out  << '\t' << (tree[it->second].groupCount[mGroups[i]]/(double)groupmap->getNumSeqs(mGroups[i]));
                            }
                        }
                    }else if (ct != NULL) {
                        if (ct->hasGroupInfo()) {
                            vector<string> mGroups = ct->getNamesOfGroups();
                            for (int i = 0; i < mGroups.size(); i++) {
                                if (((output == "detail") && (printlevel >= tree[it->second].level)) || (printlevel == tree[it->second].level)) {
                                    out  << '\t' << (tree[it->second].groupCount[mGroups[i]]/(double)ct->getGroupCount(mGroups[i]));
                                }
                            }
                        }
                    }
                }else {
                    map<string, int>::iterator itGroup;
                    if (groupmap != NULL) {
                        vector<string> mGroups = groupmap->getNamesOfGroups();
                        for (int i = 0; i < mGroups.size(); i++) {
                            if (((output == "detail") && (printlevel >= tree[it->second].level)) || (printlevel == tree[it->second].level)) {
                                out  << '\t' << tree[it->second].groupCount[mGroups[i]];
                            }
                        }
                    }else if (ct != NULL) {
                        if (ct->hasGroupInfo()) {
                            vector<string> mGroups = ct->getNamesOfGroups();
                            for (int i = 0; i < mGroups.size(); i++) {
                                if (((output == "detail") && (printlevel >= tree[it->second].level)) || (printlevel == tree[it->second].level)) {
                                    out  << '\t' << tree[it->second].groupCount[mGroups[i]];
                                }
                            }
                        }
                    }
                    
                }
                if (((output == "detail") && (printlevel >= tree[it->second].level)) || (printlevel == tree[it->second].level)) { out << endl; }
			}
			
			print(it->second, out, output);
		}
	}
	catch(exception& e) {
		m->errorOut(e, "PhyloSummary", "print");
		exit(1);
	}
}

/**************************************************************************************************/

void PhyloSummary::print(int i, ofstream& out, bool relabund){
	try {
		map<string,int>::iterator it;
		for(it=tree[i].children.begin();it!=tree[i].children.end();it++){
			
			if (tree[it->second].total != 0)  {
			
				int totalChildrenInTree = 0;
		
				map<string,int>::iterator it2;
				for(it2=tree[it->second].children.begin();it2!=tree[it->second].children.end();it2++){   
					if (tree[it2->second].total != 0)  {   totalChildrenInTree++; }
				}
                
                string nodeName = "";
                int thisNode = it->second;
                while (tree[thisNode].rank != "0") { //while you are not at top
                    if (m->getControl_pressed()) { break; }
                    nodeName = tree[thisNode].name + "|" + nodeName;
                    thisNode = tree[thisNode].parent;
                }
                if (nodeName != "") { nodeName = nodeName.substr(0, nodeName.length()-1); }
                
				out << nodeName << "\t" << (tree[it->second].total / (float)tree[i].total);
				
				map<string, int>::iterator itGroup;
				if (groupmap != NULL) {
					vector<string> mGroups = groupmap->getNamesOfGroups();
					for (int j = 0; j < mGroups.size(); j++) {
                        if (tree[i].groupCount[mGroups[j]] == 0) {
                            out  << '\t' << 0;
                        }else { out  << '\t' << (tree[it->second].groupCount[mGroups[j]] / (float)tree[i].groupCount[mGroups[j]]); }
                    }
				}else if (ct != NULL) {
                    if (ct->hasGroupInfo()) {
                        vector<string> mGroups = ct->getNamesOfGroups();
                        for (int j = 0; j < mGroups.size(); j++) {
                            if (tree[i].groupCount[mGroups[j]] == 0) {
                                out  << '\t' << 0 ;
                            }else { out  << '\t' << (tree[it->second].groupCount[mGroups[j]] / (float)tree[i].groupCount[mGroups[j]]); }
                        }
                    }
                }
				out << endl;
				
			}
			
			print(it->second, out, relabund);
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
	
		//read version
		string line = util.getline(in); util.gobble(in);
		
		int num;
		
		in >> num; util.gobble(in);
		
		tree.resize(num);
		
		in >> maxLevel; util.gobble(in);
	
		//read the tree file
		for (int i = 0; i < tree.size(); i++) {
	
			in >> tree[i].level >> num; util.gobble(in); //num contains the number of children tree[i] has
            tree[i].name = util.getline(in); util.gobble(in);
            
			//set children
			string childName;
			int childIndex;
			for (int j = 0; j < num; j++) {
				in >> childIndex; util.gobble(in);
                childName = util.getline(in); util.gobble(in);
				tree[i].children[childName] = childIndex;
			}
			
			//initialize groupcounts
			if (groupmap != NULL) {
				for (int j = 0; j < (groupmap->getNamesOfGroups()).size(); j++) {
					tree[i].groupCount[(groupmap->getNamesOfGroups())[j]] = 0;
				}
			}else if (ct != NULL) {
                if (ct->hasGroupInfo()) {
                    for (int j = 0; j < (ct->getNamesOfGroups()).size(); j++) {
                        tree[i].groupCount[(ct->getNamesOfGroups())[j]] = 0;
                    }
                }
            }
			
			tree[i].total = 0;
			
			util.gobble(in);
			
			//if (tree[i].level > maxLevel) {  maxLevel = tree[i].level;  }
		}

	}
	catch(exception& e) {
		m->errorOut(e, "PhyloSummary", "readTreeStruct");
		exit(1);
	}
}
/**************************************************************************************************/


	
