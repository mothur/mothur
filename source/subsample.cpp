//
//  subsample.cpp
//  Mothur
//
//  Created by Sarah Westcott on 4/2/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#include "subsample.h"
//**********************************************************************************************************************
Tree* SubSample::getSample(Tree* T, CountTable* ct, CountTable* newCt, int size, vector<string>& mGroups) {
    try {
        //remove seqs not in sample from counttable
        vector<string> Groups = ct->getNamesOfGroups();
        if (mGroups.size() == 0) { mGroups = Groups; }
        
        newCt->copy(ct);
        newCt->addGroup("doNotIncludeMe");
        
        map<string, int> doNotIncludeTotals; 
        vector<string> namesSeqs = ct->getNamesOfSeqs();
        for (int i = 0; i < namesSeqs.size(); i++) {  doNotIncludeTotals[namesSeqs[i]] = 0; }
    
        for (int i = 0; i < Groups.size(); i++) {
            if (util.inUsersGroups(Groups[i], mGroups)) {
                if (m->getControl_pressed()) { break; }
        
                int thisSize = ct->getGroupCount(Groups[i]);
                
                if (thisSize >= size) {	
                    
                    vector<string> names = ct->getNamesOfSeqs(Groups[i]);
                    vector<int> random;
                    for (int j = 0; j < names.size(); j++) {
                        int num = ct->getGroupCount(names[j], Groups[i]);
                        for (int k = 0; k < num; k++) { random.push_back(j); }
                    }
                    util.mothurRandomShuffle(random);
                    
                    vector<int> sampleRandoms; sampleRandoms.resize(names.size(), 0);
                    for (int j = 0; j < size; j++) { sampleRandoms[random[j]]++; }
                    for (int j = 0; j < sampleRandoms.size(); j++) {
                        newCt->setAbund(names[j], Groups[i], sampleRandoms[j]);
                    }
                    sampleRandoms.clear(); sampleRandoms.resize(names.size(), 0);
                    for (int j = size; j < thisSize; j++) { sampleRandoms[random[j]]++; }
                    for (int j = 0; j < sampleRandoms.size(); j++) {  doNotIncludeTotals[names[j]] += sampleRandoms[j]; }
                }else {  m->mothurOut("[ERROR]: You have selected a size that is larger than "+Groups[i]+" number of sequences.\n"); m->setControl_pressed(true); }
            }

        }
        
        for (map<string, int>::iterator it = doNotIncludeTotals.begin(); it != doNotIncludeTotals.end(); it++) {  
            newCt->setAbund(it->first, "doNotIncludeMe", it->second);
        } 
        
        vector<string> Treenames = T->getTreeNames();
        Tree* newTree = new Tree(newCt, Treenames);
        newTree->getCopy(T, true);
        
        return newTree;
    }
    catch(exception& e) {
        m->errorOut(e, "SubSample", "getSample-Tree");
        exit(1);
    }
}
//**********************************************************************************************************************
Tree* SubSample::getSampleWithReplacement(Tree* T, CountTable* ct, CountTable* newCt, int size, vector<string>& mGroups) {
    try {
        Tree* newTree = NULL;
        
        //remove seqs not in sample from counttable
        vector<string> Groups = ct->getNamesOfGroups();
        if (mGroups.size() == 0) { mGroups = Groups; }
        
        newCt->copy(ct);
        newCt->addGroup("doNotIncludeMe");
        
        map<string, int> doNotIncludeTotals;
        vector<string> namesSeqs = ct->getNamesOfSeqs();
        for (int i = 0; i < namesSeqs.size(); i++) {  doNotIncludeTotals[namesSeqs[i]] = 0; }
        
        for (int i = 0; i < Groups.size(); i++) {
            if (util.inUsersGroups(Groups[i], mGroups)) {
                if (m->getControl_pressed()) { break; }
                
                int thisSize = ct->getGroupCount(Groups[i]);
                
                vector<string> names = ct->getNamesOfSeqs(Groups[i]);
                vector<int> random;
                for (int j = 0; j < names.size(); j++) {
                    int num = ct->getGroupCount(names[j], Groups[i]);
                    for (int k = 0; k < num; k++) { random.push_back(j); }
                }
                
                vector<int> sampleRandoms; sampleRandoms.resize(names.size(), 0);
                long long totalNumReads = random.size()-1;
                set<long long> selected;
                for (int j = 0; j < size; j++) { //allows for multiple selection of same read
                    //"grab random from bag"
                    long long randomRead = util.getRandomIndex(totalNumReads);
                    sampleRandoms[random[randomRead]]++;
                    selected.insert(randomRead);
                }
                for (int j = 0; j < sampleRandoms.size(); j++) { //create new count file with updated sequence counts
                    newCt->setAbund(names[j], Groups[i], sampleRandoms[j]);
                }
                
                //set unselected reads to "do not include"
                sampleRandoms.clear(); sampleRandoms.resize(names.size(), 0);
                for (long long j = 0; j < random.size(); j++) {
                    if (selected.count(j) == 0)  { //we did not selected this read from random
                        sampleRandoms[random[j]]++;
                    }
                }
                for (int j = 0; j < sampleRandoms.size(); j++) {  doNotIncludeTotals[names[j]] += sampleRandoms[j]; }
            }
            
        }
        
        for (map<string, int>::iterator it = doNotIncludeTotals.begin(); it != doNotIncludeTotals.end(); it++) {
            newCt->setAbund(it->first, "doNotIncludeMe", it->second);
        }
        
        vector<string> Treenames = T->getTreeNames();
        newTree = new Tree(newCt, Treenames);
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
            
            if (m->getControl_pressed()) { break; }
            
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
            }else { m->mothurOut("[ERROR]: "+dupName+" is not in your name file, please correct.\n"); m->setControl_pressed(true); }
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

set<long long> SubSample::getWeightedSample(map<long long, long long> & nameMap, long long num) {
    try {
        set<long long> sampleNames;
    
        long long totalSeqs = nameMap.size();
        if (totalSeqs < num) { m->mothurOut("[ERROR]: Requesting sample size larger than number of seqeunces, quitting.\n"); m->setControl_pressed(true); return sampleNames; }
        else if (totalSeqs == num) {
            for (map<long long, long long>::iterator it = nameMap.begin(); it != nameMap.end(); it++) { sampleNames.insert(it->first);     }
            return sampleNames;
        }
        
        long long numSampled = 0;
        map<long long, set<long long> > weights;//weight -> names of seqs with that weight
        map<long long, set<long long> >::iterator itWeight;
        long long total = 0;
        for (map<long long, long long>::iterator it = nameMap.begin(); it != nameMap.end(); it++) {
            total += it->second;
            itWeight = weights.find(it->second);
            if (itWeight == weights.end()) { //this is a weight we haven't seen before
                set<long long> temp;
                temp.insert(it->first);
                weights[it->second] = temp;
            }else {
                weights[it->second].insert(it->first); //dup weight, combine to save memory
            }
        }
        
        //find running total
        long long runningTotal = 0;
        map<long long, set<long long> > cumulative;//weight + sum so far -> names of seqs with that weight
        for (itWeight = weights.begin(); itWeight != weights.end(); itWeight++) {
            int count = itWeight->second.size(); //number of seqs with this weight
            runningTotal += itWeight->first * count;
            cumulative[runningTotal] = itWeight->second;
        }
        weights.clear();
        
        while(numSampled != num) {
            long long index = util.getRandomIndex(total); //random index including weights
            
            map<long long, set<long long> >::iterator itWeight = cumulative.lower_bound(index);
            
            sampleNames.insert(*itWeight->second.begin()); //save name in sample names
            
            itWeight->second.erase(itWeight->second.begin()); //remove seq since we sampled it
            
            if (itWeight->second.size() == 0) { cumulative.erase(itWeight); total = cumulative.rbegin()->first;  } //remove this weight if all seqs are sampled. Reset bound if needed.
            
            numSampled = sampleNames.size();
        }
    
        return sampleNames;
    }
    catch(exception& e) {
        m->errorOut(e, "SubSample", "getWeightedSample");
        exit(1);
    }
}
//**********************************************************************************************************************
vector<string> SubSample::getSample(vector<SharedRAbundVector*>& rabunds, int size, vector<string> currentLabels) {
    try {
        
        //save mothurOut's binLabels to restore for next label
        vector<string> saveBinLabels = currentLabels;
        SharedRAbundVectors* newLookup = new SharedRAbundVectors();
        
        int numBins = rabunds[0]->getNumBins();
        for (int i = 0; i < rabunds.size(); i++) {
            int thisSize = rabunds[i]->getNumSeqs();
            
            if (thisSize != size) {
                
                vector<int> order;
                for (int j = 0; j < rabunds[i]->size(); j++) {
                    int abund = rabunds[i]->get(j);
                    for(int k=0;k<abund;k++){ order.push_back(j);  }
                }
                
                util.mothurRandomShuffle(order);
                
                SharedRAbundVector* temp = new SharedRAbundVector(numBins);
                temp->setLabel(rabunds[i]->getLabel());
                temp->setGroup(rabunds[i]->getGroup());
                
                for (int j = 0; j < size; j++) { //only allows you to select a read once
                    if (m->getControl_pressed()) {  return currentLabels; }
                    temp->increment(order[j]);
                }
                newLookup->push_back(temp);
            }else { SharedRAbundVector* temp = new SharedRAbundVector(*rabunds[i]); newLookup->push_back(temp); }
        }
        newLookup->setOTUNames(currentLabels);
        newLookup->eliminateZeroOTUS();
        
        for (int i = 0; i < rabunds.size(); i++) { delete rabunds[i]; } rabunds.clear();
        rabunds = newLookup->getSharedRAbundVectors();
       
        //save mothurOut's binLabels to restore for next label
        vector<string> subsampleBinLabels = newLookup->getOTUNames();
        delete newLookup;
        
        return subsampleBinLabels;
        
    }
    catch(exception& e) {
        m->errorOut(e, "SubSample", "getSample-shared");
        exit(1);
    }
}
//**********************************************************************************************************************
vector<string> SubSample::getSampleWithReplacement(vector<SharedRAbundVector*>& rabunds, int size, vector<string> currentLabels) {
    try {
        
        //save mothurOut's binLabels to restore for next label
        vector<string> saveBinLabels = currentLabels;
        SharedRAbundVectors* newLookup = new SharedRAbundVectors();
        
        int numBins = rabunds[0]->getNumBins();
        for (int i = 0; i < rabunds.size(); i++) {
            int thisSize = rabunds[i]->getNumSeqs();

            vector<int> order;
            for (int j = 0; j < rabunds[i]->size(); j++) {
                int abund = rabunds[i]->get(j);
                for(int k=0;k<abund;k++){ order.push_back(j);  }
            }
            
            SharedRAbundVector* temp = new SharedRAbundVector(numBins);
            temp->setLabel(rabunds[i]->getLabel());
            temp->setGroup(rabunds[i]->getGroup());
            
            long long orderSize = order.size()-1;
            for (int j = 0; j < size; j++) { //allows you to select a read multiple times
                if (m->getControl_pressed()) {  return currentLabels; }
                //"grab random from bag"
                long long randomRead = util.getRandomIndex(orderSize);
                temp->increment(order[randomRead]);
            }
            newLookup->push_back(temp);

        }
        newLookup->setOTUNames(currentLabels);
        newLookup->eliminateZeroOTUS();
        
        for (int i = 0; i < rabunds.size(); i++) { delete rabunds[i]; } rabunds.clear();
        rabunds = newLookup->getSharedRAbundVectors();
        
        //save mothurOut's binLabels to restore for next label
        vector<string> subsampleBinLabels = newLookup->getOTUNames();
        delete newLookup;
        
        return subsampleBinLabels;
        
    }
    catch(exception& e) {
        m->errorOut(e, "SubSample", "getSampleWithReplacement-shared");
        exit(1);
    }
}
//**********************************************************************************************************************
vector<string> SubSample::getSample(SharedRAbundVectors*& thislookup, int size) {
	try {
		//save mothurOut's binLabels to restore for next label
		vector<string> saveBinLabels = thislookup->getOTUNames();
        vector<SharedRAbundVector*> rabunds = thislookup->getSharedRAbundVectors();
        
        vector<string> subsampleBinLabels = getSample(rabunds, size, saveBinLabels);
        SharedRAbundVectors* newLookup = new SharedRAbundVectors();
		
		for (int i = 0; i < rabunds.size(); i++) { newLookup->push_back(rabunds[i]);  }
        newLookup->setOTUNames(subsampleBinLabels);
        delete thislookup;
        thislookup = newLookup;
        
		return subsampleBinLabels;
		
	}
	catch(exception& e) {
		m->errorOut(e, "SubSample", "getSample-shared");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> SubSample::getSampleWithReplacement(SharedRAbundVectors*& thislookup, int size) {
    try {
        //save mothurOut's binLabels to restore for next label
        vector<string> saveBinLabels = thislookup->getOTUNames();
        vector<SharedRAbundVector*> rabunds = thislookup->getSharedRAbundVectors();
        
        vector<string> subsampleBinLabels = getSampleWithReplacement(rabunds, size, saveBinLabels);
        SharedRAbundVectors* newLookup = new SharedRAbundVectors();
        
        for (int i = 0; i < rabunds.size(); i++) { newLookup->push_back(rabunds[i]);  }
        newLookup->setOTUNames(subsampleBinLabels);
        delete thislookup;
        thislookup = newLookup;
        
        return subsampleBinLabels;
        
    }
    catch(exception& e) {
        m->errorOut(e, "SubSample", "getSampleWithReplacement-shared");
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
			util.mothurRandomShuffle(order);
			
            RAbundVector rabund(numBins);
			rabund.setLabel(sabund->getLabel());

			for (int j = 0; j < size; j++) {
                
				if (m->getControl_pressed()) { return 0; }
				
                int abund = rabund.get(order.get(j));
				rabund.set(order.get(j), (abund+1));
			}

            delete sabund;
            sabund = new SAbundVector();
            *sabund = rabund.getSAbundVector();
            
		}else if (thisSize < size) { m->mothurOut("[ERROR]: The size you requested is larger than the number of sequences in the sabund vector. You requested " + toString(size) + " and you only have " + toString(thisSize) + " seqs in your sabund vector.\n"); m->setControl_pressed(true); }
        
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "SubSampleCommand", "getSample");
		exit(1);
	}
}
//**********************************************************************************************************************
int SubSample::getSampleWithReplacement(SAbundVector*& sabund, int size) {
    try {
        
        int numBins = sabund->getNumBins();
        int thisSize = sabund->getNumSeqs();
        
        OrderVector order = sabund->getOrderVector();
        
        if (thisSize > size) {
            RAbundVector rabund(numBins);
            rabund.setLabel(sabund->getLabel());
            
            long long orderSize = order.size()-1;
            for (int j = 0; j < size; j++) {
                
                if (m->getControl_pressed()) { return 0; }
                
                //"grab random from bag"
                long long randomRead = util.getRandomIndex(orderSize);
                
                int abund = rabund.get(order.get(randomRead));
                rabund.set(order.get(randomRead), (abund+1));
            }
            
            delete sabund;
            sabund = new SAbundVector();
            *sabund = rabund.getSAbundVector();
            
        }else if (thisSize < size) { m->mothurOut("[ERROR]: The size you requested is larger than the number of sequences in the sabund vector. You requested " + toString(size) + " and you only have " + toString(thisSize) + " seqs in your sabund vector.\n"); m->setControl_pressed(true); }
        
        return 0;
        
    }
    catch(exception& e) {
        m->errorOut(e, "SubSampleCommand", "getSampleWithReplacement");
        exit(1);
    }
}
//**********************************************************************************************************************
int SubSample::getSample(RAbundVector*& rabund, int size) {
    try {
        
        int numBins = rabund->getNumBins();
        int thisSize = rabund->getNumSeqs();
        
        OrderVector order = rabund->getOrderVector(NULL);
        
        if (thisSize > size) {
            util.mothurRandomShuffle(order);
            
            RAbundVector sampledRabund(numBins);
            sampledRabund.setLabel(rabund->getLabel());
            
            for (int j = 0; j < size; j++) {
                
                if (m->getControl_pressed()) { return 0; }
                
                int abund = sampledRabund.get(order.get(j));
                sampledRabund.set(order.get(j), (abund+1));
            }
            
            delete rabund;
            rabund = new RAbundVector(sampledRabund);
            
        }else if (thisSize < size) { m->mothurOut("[ERROR]: The size you requested is larger than the number of sequences in the rabund vector. You requested " + toString(size) + " and you only have " + toString(thisSize) + " seqs in your rabund vector.\n"); m->setControl_pressed(true); }
        
        return 0;
        
    }
    catch(exception& e) {
        m->errorOut(e, "SubSampleCommand", "getSample");
        exit(1);
    }
}
//**********************************************************************************************************************
int SubSample::getSampleWithReplacement(RAbundVector*& rabund, int size) {
    try {
        
        int numBins = rabund->getNumBins();
        int thisSize = rabund->getNumSeqs();
        
        OrderVector order = rabund->getOrderVector(NULL);
        
        if (thisSize > size) {
            RAbundVector sampledRabund(numBins);
            sampledRabund.setLabel(rabund->getLabel());
            
            long long orderSize = order.size()-1;
            for (int j = 0; j < size; j++) {
                
                if (m->getControl_pressed()) { return 0; }
                
                //"grab random from bag"
                long long randomRead = util.getRandomIndex(orderSize);
                
                int abund = sampledRabund.get(order.get(randomRead));
                sampledRabund.set(order.get(randomRead), (abund+1));
            }
            
            delete rabund;
            rabund = new RAbundVector(sampledRabund);
            
        }else if (thisSize < size) { m->mothurOut("[ERROR]: The size you requested is larger than the number of sequences in the sabund vector. You requested " + toString(size) + " and you only have " + toString(thisSize) + " seqs in your sabund vector.\n"); m->setControl_pressed(true); }
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "SubSampleCommand", "getSampleWithReplacement");
        exit(1);
    }
}
//**********************************************************************************************************************
CountTable SubSample::getSample(CountTable& ct, int size, vector<string> Groups, bool persample) {
	try {
        if (!persample) { return (getSample(ct, size, Groups)); }
        
        if (!ct.hasGroupInfo()) { m->mothurOut("[ERROR]: Cannot subsample by group because your count table doesn't have group information.\n"); m->setControl_pressed(true); }
            
        CountTable sampledCt;
        map<string, vector<int> > tempCount;
        for (int i = 0; i < Groups.size(); i++) {
            sampledCt.addGroup(Groups[i]);
            
            vector<string> names = ct.getNamesOfSeqs(Groups[i]);
            vector<int> allNames;
            for (int j = 0; j < names.size(); j++) {
                
                if (m->getControl_pressed()) { return sampledCt; }
                
                int num = ct. getGroupCount(names[j], Groups[i]);
                for (int k = 0; k < num; k++) { allNames.push_back(j); }
            }
            
            util.mothurRandomShuffle(allNames);
            
            if (allNames.size() < size) { m->mothurOut("[ERROR]: You have selected a size that is larger than "+Groups[i]+" number of sequences.\n"); m->setControl_pressed(true); }
            else{
                for (int j = 0; j < size; j++) {
                    
                    if (m->getControl_pressed()) { return sampledCt; }
                    
                    map<string, vector<int> >::iterator it = tempCount.find(names[allNames[j]]);
                    
                    if (it == tempCount.end()) { //we have not seen this sequence at all yet
                        vector<int> tempGroups; tempGroups.resize(Groups.size(), 0);
                        tempGroups[i]++;
                        tempCount[names[allNames[j]]] = tempGroups;
                    }else{
                        tempCount[names[allNames[j]]][i]++;
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
CountTable SubSample::getSampleWithReplacement(CountTable& ct, int size, vector<string> Groups, bool persample) {
    try {
        if (!persample) { return (getSampleWithReplacement(ct, size, Groups)); }
        
        if (!ct.hasGroupInfo()) { m->mothurOut("[ERROR]: Cannot subsample by group because your count table doesn't have group information.\n"); m->setControl_pressed(true); }
        
        CountTable sampledCt;
        map<string, vector<int> > tempCount;
        for (int i = 0; i < Groups.size(); i++) {
            sampledCt.addGroup(Groups[i]);
            
            vector<string> names = ct.getNamesOfSeqs(Groups[i]);
            vector<int> allNames;
            for (int j = 0; j < names.size(); j++) {
                
                if (m->getControl_pressed()) { return sampledCt; }
                
                int num = ct.getGroupCount(names[j], Groups[i]);
                for (int k = 0; k < num; k++) { allNames.push_back(j); }
            }
            
            long long allNamesSize = allNames.size()-1;
            
            if (allNames.size() < size) { m->mothurOut("[ERROR]: You have selected a size that is larger than "+Groups[i]+" number of sequences.\n"); m->setControl_pressed(true); }
            else{
                for (int j = 0; j < size; j++) {
                    
                    if (m->getControl_pressed()) { return sampledCt; }
                    
                    //"grab random from bag"
                    long long randomRead = util.getRandomIndex(allNamesSize);
                    
                    map<string, vector<int> >::iterator it = tempCount.find(names[allNames[randomRead]]);
                    
                    if (it == tempCount.end()) { //we have not seen this sequence at all yet
                        vector<int> tempGroups; tempGroups.resize(Groups.size(), 0);
                        tempGroups[i]++;
                        tempCount[names[allNames[randomRead]]] = tempGroups;
                    }else{
                        tempCount[names[allNames[randomRead]]][i]++;
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
        m->errorOut(e, "SubSampleCommand", "getSampleWithReplacement");
        exit(1);
    }
}
//**********************************************************************************************************************
GroupMap SubSample::getSample(GroupMap& groupMap, int size, vector<string> Groups, bool persample) {
    try {
        if (!persample) { return (getSample(groupMap, size, Groups)); }
        
        GroupMap sampledGM;
        
        //initialize counts
        map<string, int> groupCounts;
        map<string, int>::iterator itGroupCounts;
        for (int i = 0; i < Groups.size(); i++) { groupCounts[Groups[i]] = 0; }
        
        for (int i = 0; i < Groups.size(); i++) {
            
            if (m->getControl_pressed()) { break; }
            
            string thisGroup = Groups[i];
            int thisSize = groupMap.getNumSeqs(thisGroup);
            
            if (thisSize >= size) {
                
                vector<string> names = groupMap.getNamesSeqs(thisGroup);
                
                util.mothurRandomShuffle(names);
                
                for (int j = 0; j < size; j++) { sampledGM.addSeq(names[j], thisGroup); }
                
            }else {  m->mothurOut("[ERROR]: You have selected a size that is larger than "+Groups[i]+" number of sequences.\n"); m->setControl_pressed(true); }
        }

        return sampledGM;
    }
    catch(exception& e) {
        m->errorOut(e, "SubSampleCommand", "getSample");
        exit(1);
    }
}
//**********************************************************************************************************************
GroupMap SubSample::getSample(GroupMap& groupMap, int size, vector<string> Groups) {
    try {
        GroupMap sampledGM;
        
        int thisSize = groupMap.getNumSeqs();
        
        if (thisSize >= size) {
            
            vector<string> names = groupMap.getNamesSeqs();
            
            util.mothurRandomShuffle(names);
            
            int numSelected = 0;
            for (int j = 0; j < names.size(); j++) {
                string thisGroup = groupMap.getGroup(names[j]);
                
                if (util.inUsersGroups(thisGroup, Groups)) {
                    sampledGM.addSeq(names[j], thisGroup);
                    numSelected++;
                }
                
                //do we have enough??
                if (numSelected == size) { break; }
            }
            
        }else {  m->mothurOut("[ERROR]: You have selected a size that is larger than the number of sequences.\n"); m->setControl_pressed(true); }
        
        return sampledGM;
    }
    catch(exception& e) {
        m->errorOut(e, "SubSampleCommand", "getSample");
        exit(1);
    }
}
//**********************************************************************************************************************
GroupMap SubSample::getSample(GroupMap& groupMap, int size) {
    try {
        GroupMap sampledGM;
    
        int thisSize = groupMap.getNumSeqs();
        
        if (thisSize >= size) {
            
            vector<string> names = groupMap.getNamesSeqs();
            
            util.mothurRandomShuffle(names);
            
            for (int j = 0; j < size; j++) {
                string thisGroup = groupMap.getGroup(names[j]);
                sampledGM.addSeq(names[j], thisGroup);
            }
            
        }else {  m->mothurOut("[ERROR]: You have selected a size that is larger than the number of sequences.\n"); m->setControl_pressed(true); }
        
        return sampledGM;
    }
    catch(exception& e) {
        m->errorOut(e, "SubSampleCommand", "getSample");
        exit(1);
    }
}
//**********************************************************************************************************************
CountTable SubSample::getSample(CountTable& ct, int size, vector<string> Groups) {
	try {
        CountTable sampledCt;
        if (!ct.hasGroupInfo() && (Groups.size() != 0)) { m->mothurOut("[ERROR]: Cannot subsample with groups because your count table doesn't have group information.\n"); m->setControl_pressed(true); return sampledCt; }
        
        if (ct.hasGroupInfo()) { //only select reads from Groups
        
            map<string, vector<int> > tempCount;
            for (int i = 0; i < Groups.size(); i++) { sampledCt.addGroup(Groups[i]);  }
                
            vector<string> names = ct.getNamesOfSeqs(Groups); //names of sequences in groups
            vector<item> allNames;
            for (int j = 0; j < names.size(); j++) {
                    
                if (m->getControl_pressed()) { return sampledCt; }
                
                for (int i = 0; i < Groups.size(); i++) {
                    int num = ct.getGroupCount(names[j], Groups[i]); //num reads in this group from this seq
                    item thisSeq(j,i);
                    for (int k = 0; k < num; k++) { allNames.push_back(thisSeq); }
                }
            }
                
            util.mothurRandomShuffle(allNames);
                
            if (allNames.size() < size) { m->mothurOut("[ERROR]: You have selected a size that is larger than the number of sequences.\n"); m->setControl_pressed(true); }
            else{
                for (int j = 0; j < size; j++) {
                        
                    if (m->getControl_pressed()) { return sampledCt; }
                        
                    map<string, vector<int> >::iterator it = tempCount.find(names[allNames[j].name]);
                        
                    if (it == tempCount.end()) { //we have not seen this sequence at all yet
                        vector<int> tempGroups; tempGroups.resize(Groups.size(), 0);
                        tempGroups[allNames[j].group]++;
                        tempCount[names[allNames[j].name]] = tempGroups;
                    }else{
                        tempCount[names[allNames[j].name]][allNames[j].group]++;
                    }
                }
            }
            
            
            //build count table
            for (map<string, vector<int> >::iterator it = tempCount.begin(); it != tempCount.end();) {
                sampledCt.push_back(it->first, it->second);
                tempCount.erase(it++);
            }

        }else { //no groups
            vector<string> names = ct.getNamesOfSeqs();
            map<int, int> nameMap;
            vector<int> allNames;
            
            for (int i = 0; i < names.size(); i++) {
                int num = ct.getNumSeqs(names[i]);
                for (int j = 0; j < num; j++) { allNames.push_back(i); }
            }
            
            if (allNames.size() < size) { m->mothurOut("[ERROR]: You have selected a size that is larger than the number of sequences.\n"); m->setControl_pressed(true); return sampledCt; }
            else {
                util.mothurRandomShuffle(allNames);
                
                for (int j = 0; j < size; j++) {
                    if (m->getControl_pressed()) { return sampledCt; }
                    
                    map<int, int>::iterator it = nameMap.find(allNames[j]);
                    
                    //we have not seen this sequence at all yet
                    if (it == nameMap.end()) { nameMap[allNames[j]] = 1;  }
                    else{  nameMap[allNames[j]]++;  }
                }
                
                //build count table
                for (map<int, int>::iterator it = nameMap.begin(); it != nameMap.end();) {
                    sampledCt.push_back(names[it->first], it->second);
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
CountTable SubSample::getSampleWithReplacement(CountTable& ct, int size, vector<string> Groups) {
    try {
        CountTable sampledCt;
        if (!ct.hasGroupInfo() && (Groups.size() != 0)) { m->mothurOut("[ERROR]: Cannot subsample with groups because your count table doesn't have group information.\n"); m->setControl_pressed(true); return sampledCt; }
        
        if (ct.hasGroupInfo()) { //only select reads from Groups
            
            map<string, vector<int> > tempCount;
            for (int i = 0; i < Groups.size(); i++) { sampledCt.addGroup(Groups[i]);  }
            
            vector<string> names = ct.getNamesOfSeqs(Groups); //names of sequences in groups
            vector<item> allNames;
            for (int j = 0; j < names.size(); j++) {
                
                if (m->getControl_pressed()) { return sampledCt; }
                
                for (int i = 0; i < Groups.size(); i++) {
                    int num = ct.getGroupCount(names[j], Groups[i]); //num reads in this group from this seq
                    item thisSeq(j,i);
                    for (int k = 0; k < num; k++) { allNames.push_back(thisSeq); }
                }
            }
            
            if (allNames.size() < size) { m->mothurOut("[ERROR]: You have selected a size that is larger than the number of sequences.\n"); m->setControl_pressed(true); }
            else{
                long long allNamesSize = allNames.size()-1;
                
                for (int j = 0; j < size; j++) {
                    
                    if (m->getControl_pressed()) { return sampledCt; }
                    
                    long long randomRead = util.getRandomIndex(allNamesSize);
                    
                    map<string, vector<int> >::iterator it = tempCount.find(names[allNames[randomRead].name]);
                    
                    if (it == tempCount.end()) { //we have not seen this sequence at all yet
                        vector<int> tempGroups; tempGroups.resize(Groups.size(), 0);
                        tempGroups[allNames[randomRead].group]++;
                        tempCount[names[allNames[randomRead].name]] = tempGroups;
                    }else{
                        tempCount[names[allNames[randomRead].name]][allNames[randomRead].group]++;
                    }
                }
            }
            
            
            //build count table
            for (map<string, vector<int> >::iterator it = tempCount.begin(); it != tempCount.end();) {
                sampledCt.push_back(it->first, it->second);
                tempCount.erase(it++);
            }
            
        }else { //no groups
            vector<string> names = ct.getNamesOfSeqs();
            map<int, int> nameMap;
            vector<int> allNames;
            
            for (int i = 0; i < names.size(); i++) {
                int num = ct.getNumSeqs(names[i]);
                for (int j = 0; j < num; j++) { allNames.push_back(i); }
            }
            
            if (allNames.size() < size) { m->mothurOut("[ERROR]: You have selected a size that is larger than the number of sequences.\n"); m->setControl_pressed(true); return sampledCt; }
            else {
                long long allNamesSize = allNames.size()-1;
                
                for (int j = 0; j < size; j++) {
                    if (m->getControl_pressed()) { return sampledCt; }
                    
                    long long randomRead = util.getRandomIndex(allNamesSize);
                    
                    map<int, int>::iterator it = nameMap.find(allNames[randomRead]);
                    
                    //we have not seen this sequence at all yet
                    if (it == nameMap.end()) { nameMap[allNames[randomRead]] = 1;  }
                    else{  nameMap[allNames[randomRead]]++;  }
                }
                
                //build count table
                for (map<int, int>::iterator it = nameMap.begin(); it != nameMap.end();) {
                    sampledCt.push_back(names[it->first], it->second);
                    nameMap.erase(it++);
                }
            }
        }
        
        return sampledCt;
    }
    catch(exception& e) {
        m->errorOut(e, "SubSampleCommand", "getSampleWithReplacement");
        exit(1);
    }
}
//**********************************************************************************************************************


