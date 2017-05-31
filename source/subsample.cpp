//
//  subsample.cpp
//  Mothur
//
//  Created by Sarah Westcott on 4/2/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#include "subsample.h"
//**********************************************************************************************************************
Tree* SubSample::getSample(Tree* T, CountTable* ct, CountTable* newCt, int size) {
    try {
        Tree* newTree = NULL;
        
        //remove seqs not in sample from counttable
        vector<string> Groups = ct->getNamesOfGroups();
        newCt->copy(ct);
        newCt->addGroup("doNotIncludeMe");
        
        map<string, int> doNotIncludeTotals; 
        vector<string> namesSeqs = ct->getNamesOfSeqs();
        for (int i = 0; i < namesSeqs.size(); i++) {  doNotIncludeTotals[namesSeqs[i]] = 0; }
    
        for (int i = 0; i < Groups.size(); i++) {
            if (m->inUsersGroups(Groups[i], m->getGroups())) {
                if (m->control_pressed) { break; }
        
                int thisSize = ct->getGroupCount(Groups[i]);
                
                if (thisSize >= size) {	
                    
                    vector<string> names = ct->getNamesOfSeqs(Groups[i]);
                    vector<int> random;
                    for (int j = 0; j < names.size(); j++) {
                        int num = ct->getGroupCount(names[j], Groups[i]);
                        for (int k = 0; k < num; k++) { random.push_back(j); }
                    }
                    m->mothurRandomShuffle(random);
                    
                    vector<int> sampleRandoms; sampleRandoms.resize(names.size(), 0);
                    for (int j = 0; j < size; j++) { sampleRandoms[random[j]]++; }
                    for (int j = 0; j < sampleRandoms.size(); j++) {
                        newCt->setAbund(names[j], Groups[i], sampleRandoms[j]);
                    }
                    sampleRandoms.clear(); sampleRandoms.resize(names.size(), 0);
                    for (int j = size; j < thisSize; j++) { sampleRandoms[random[j]]++; }
                    for (int j = 0; j < sampleRandoms.size(); j++) {  doNotIncludeTotals[names[j]] += sampleRandoms[j]; }
                }else {  m->mothurOut("[ERROR]: You have selected a size that is larger than "+Groups[i]+" number of sequences.\n"); m->control_pressed = true; }
            }

        }
        
        for (map<string, int>::iterator it = doNotIncludeTotals.begin(); it != doNotIncludeTotals.end(); it++) {  
            newCt->setAbund(it->first, "doNotIncludeMe", it->second);
        } 
        
       
        newTree = new Tree(newCt);
        newTree->getCopy(T, true);
        
        return newTree;
    }
    catch(exception& e) {
        m->errorOut(e, "SubSample", "getSample-Tree");
        exit(1);
    }
}
//**********************************************************************************************************************
//assumes whole maps dupName -> uniqueName
map<string, string> SubSample::deconvolute(map<string, string> whole, vector<string>& wanted) {
    try {
        map<string, string> nameMap;
        
        //whole will be empty if user gave no name file, so we don't need to make a new one
        if (whole.size() == 0) { return nameMap; }
        
        vector<string> newWanted;
        for (int i = 0; i < wanted.size(); i++) {
            
            if (m->control_pressed) { break; }
            
            string dupName = wanted[i];
            
            map<string, string>::iterator itWhole = whole.find(dupName);
            if (itWhole != whole.end()) {
                string repName = itWhole->second;
                
                //do we already have this rep?
                map<string, string>::iterator itName = nameMap.find(repName);
                if (itName != nameMap.end()) { //add this seqs to dups list
                    (itName->second) += "," + dupName;
                }else { //first sighting of this seq
                    nameMap[repName] = dupName;
                    newWanted.push_back(repName);
                }
            }else { m->mothurOut("[ERROR]: "+dupName+" is not in your name file, please correct.\n"); m->control_pressed = true; }
        }
        
        wanted = newWanted;
        return nameMap;
    }
	catch(exception& e) {
		m->errorOut(e, "SubSample", "deconvolute");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> SubSample::getSample(SharedRAbundVectors* thislookup, int size) {
	try {
		
		//save mothurOut's binLabels to restore for next label
		vector<string> saveBinLabels = m->currentSharedBinLabels;
        vector<string> groups = thislookup->getNamesGroups();
        vector<RAbundVector*> rabunds = thislookup->getSharedRAbundVectors();
        SharedRAbundVectors* newLookup = new SharedRAbundVectors();
		
		int numBins = thislookup->getNumBins();
		for (int i = 0; i < rabunds.size(); i++) {
			int thisSize = rabunds[i]->getNumSeqs();
			
			if (thisSize != size) {
				
				OrderVector order = rabunds[i]->getOrderVector(NULL);
                
                m->mothurRandomShuffle(order);
				
				RAbundVector* temp = new RAbundVector(numBins);
				temp->setLabel(rabunds[i]->getLabel());
				
				for (int j = 0; j < size; j++) {
					
					if (m->control_pressed) {  return m->currentSharedBinLabels; }
					
					int bin = order.get(j);
					
					int abund = temp->get(bin);
					temp->set(bin, (abund+1));
				}
                newLookup->push_back(temp, groups[i]);
			}
		}
        newLookup->eliminateZeroOTUS();
        delete thislookup;
        thislookup = newLookup;
        
		if (m->control_pressed) { return m->currentSharedBinLabels; }
		
		//save mothurOut's binLabels to restore for next label
        vector<string> subsampleBinLabels = m->currentSharedBinLabels;
		m->currentSharedBinLabels = saveBinLabels;
		
		return subsampleBinLabels;
		
	}
	catch(exception& e) {
		m->errorOut(e, "SubSample", "getSample-shared");
		exit(1);
	}
}	
//**********************************************************************************************************************
int SubSample::getSample(SAbundVector*& sabund, int size) {
	try {
		
        int numBins = sabund->getNumBins();
        int thisSize = sabund->getNumSeqs();

        OrderVector order = sabund->getOrderVector();
        
		if (thisSize > size) {
			m->mothurRandomShuffle(order);
			
            RAbundVector rabund(numBins);
			rabund.setLabel(sabund->getLabel());

			for (int j = 0; j < size; j++) {
                
				if (m->control_pressed) { return 0; }
				
                int abund = rabund.get(order.get(j));
				rabund.set(order.get(j), (abund+1));
			}

            delete sabund;
            sabund = new SAbundVector();
            *sabund = rabund.getSAbundVector();
            
		}else if (thisSize < size) { m->mothurOut("[ERROR]: The size you requested is larger than the number of sequences in the sabund vector. You requested " + toString(size) + " and you only have " + toString(thisSize) + " seqs in your sabund vector.\n"); m->control_pressed = true; }
		
		if (m->control_pressed) { return 0; }
        
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "SubSampleCommand", "getSample");
		exit(1);
	}
}
//**********************************************************************************************************************
CountTable SubSample::getSample(CountTable& ct, int size, vector<string> Groups) {
	try {
        if (!ct.hasGroupInfo()) { m->mothurOut("[ERROR]: Cannot subsample by group because your count table doesn't have group information.\n"); m->control_pressed = true; }
            
        CountTable sampledCt;
        map<string, vector<int> > tempCount;
        for (int i = 0; i < Groups.size(); i++) {
            sampledCt.addGroup(Groups[i]);
            
            vector<string> names = ct.getNamesOfSeqs(Groups[i]);
            vector<string> allNames;
            for (int j = 0; j < names.size(); j++) {
                
                if (m->control_pressed) { return sampledCt; }
                
                int num = ct. getGroupCount(names[j], Groups[i]);
                for (int k = 0; k < num; k++) { allNames.push_back(names[j]); }
            }
            
            m->mothurRandomShuffle(allNames);
            
            if (allNames.size() < size) { m->mothurOut("[ERROR]: You have selected a size that is larger than "+Groups[i]+" number of sequences.\n"); m->control_pressed = true; }
            else{
                for (int j = 0; j < size; j++) {
                    
                    if (m->control_pressed) { return sampledCt; }
                    
                    map<string, vector<int> >::iterator it = tempCount.find(allNames[j]);
                    
                    if (it == tempCount.end()) { //we have not seen this sequence at all yet
                        vector<int> tempGroups; tempGroups.resize(Groups.size(), 0);
                        tempGroups[i]++;
                        tempCount[allNames[j]] = tempGroups;
                    }else{
                        tempCount[allNames[j]][i]++;
                    }
                }
            }
        }
        
        //build count table
        for (map<string, vector<int> >::iterator it = tempCount.begin(); it != tempCount.end();) {
            sampledCt.push_back(it->first, it->second);
            tempCount.erase(it++);
        }
        
        return sampledCt;
    }
	catch(exception& e) {
		m->errorOut(e, "SubSampleCommand", "getSample");
		exit(1);
	}
}
//**********************************************************************************************************************
CountTable SubSample::getSample(CountTable& ct, int size, vector<string> Groups, bool pickedGroups) {
	try {
        CountTable sampledCt;
        if (!ct.hasGroupInfo() && pickedGroups) { m->mothurOut("[ERROR]: Cannot subsample with groups because your count table doesn't have group information.\n"); m->control_pressed = true; return sampledCt; }
        
        if (ct.hasGroupInfo()) {
            map<string, vector<int> > tempCount;
            vector<item> allNames;
            map<string, int> groupMap;
            
            vector<string> myGroups;
            if (pickedGroups) { myGroups = Groups; }
            else {  myGroups = ct.getNamesOfGroups(); }
            
            for (int i = 0; i < myGroups.size(); i++) {
                sampledCt.addGroup(myGroups[i]);
                groupMap[myGroups[i]] = i;
                
                vector<string> names = ct.getNamesOfSeqs(myGroups[i]);
                for (int j = 0; j < names.size(); j++) {
                    
                    if (m->control_pressed) { return sampledCt; }
                    
                    int num = ct. getGroupCount(names[j], myGroups[i]);
                    for (int k = 0; k < num; k++) { 
                        item temp(names[j], myGroups[i]);
                        allNames.push_back(temp); 
                    }
                }
            }
            
            m->mothurRandomShuffle(allNames);
            
            if (allNames.size() < size) { 
                if (pickedGroups) { m->mothurOut("[ERROR]: You have selected a size that is larger than the number of sequences.\n"); } 
                else { m->mothurOut("[ERROR]: You have selected a size that is larger than the number of sequences in the groups you chose.\n"); }
                m->control_pressed = true; return sampledCt; }
            else{
                for (int j = 0; j < size; j++) {
                    
                    if (m->control_pressed) { return sampledCt; }
                    
                    map<string, vector<int> >::iterator it = tempCount.find(allNames[j].name);
                    
                    if (it == tempCount.end()) { //we have not seen this sequence at all yet
                        vector<int> tempGroups; tempGroups.resize(myGroups.size(), 0);
                        tempGroups[groupMap[allNames[j].group]]++;
                        tempCount[allNames[j].name] = tempGroups;
                    }else{
                        tempCount[allNames[j].name][groupMap[allNames[j].group]]++;
                    }
                }
            }
            
            //build count table
            for (map<string, vector<int> >::iterator it = tempCount.begin(); it != tempCount.end();) {
                sampledCt.push_back(it->first, it->second);
                tempCount.erase(it++);
            }
            
            //remove empty groups 
            for (int i = 0; i < myGroups.size(); i++) { if (sampledCt.getGroupCount(myGroups[i]) == 0) { sampledCt.removeGroup(myGroups[i]); } }
            
        }else {
            vector<string> names = ct.getNamesOfSeqs();
            map<string, int> nameMap;
            vector<string> allNames;
            
            for (int i = 0; i < names.size(); i++) {
                int num = ct.getNumSeqs(names[i]);
                for (int j = 0; j < num; j++) { allNames.push_back(names[i]); }
            }
            
            if (allNames.size() < size) { m->mothurOut("[ERROR]: You have selected a size that is larger than the number of sequences.\n"); m->control_pressed = true; return sampledCt; }
            else {
                m->mothurRandomShuffle(allNames);
                
                for (int j = 0; j < size; j++) {
                    if (m->control_pressed) { return sampledCt; }
                    
                    map<string, int>::iterator it = nameMap.find(allNames[j]);
                    
                    //we have not seen this sequence at all yet
                    if (it == nameMap.end()) { nameMap[allNames[j]] = 1;  }
                    else{  nameMap[allNames[j]]++;  }
                }
                
                //build count table
                for (map<string, int>::iterator it = nameMap.begin(); it != nameMap.end();) {
                    sampledCt.push_back(it->first, it->second);
                    nameMap.erase(it++);
                }
            }
        }
        
        return sampledCt;
    }
	catch(exception& e) {
		m->errorOut(e, "SubSampleCommand", "getSample");
		exit(1);
	}
}
//**********************************************************************************************************************


