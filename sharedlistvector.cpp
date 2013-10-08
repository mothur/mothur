/*
 *  sharedSharedListVector.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 1/22/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "sabundvector.hpp"
#include "rabundvector.hpp"
#include "ordervector.hpp"
#include "sharedlistvector.h"
#include "sharedordervector.h"
#include "sharedutilities.h"

/***********************************************************************/

SharedListVector::SharedListVector() : DataVector(), maxRank(0), numBins(0), numSeqs(0){ groupmap = NULL; countTable = NULL; }

/***********************************************************************/

SharedListVector::SharedListVector(int n):	DataVector(), data(n, "") , maxRank(0), numBins(0), numSeqs(0){ groupmap = NULL; countTable = NULL; }

/***********************************************************************/
SharedListVector::SharedListVector(ifstream& f) : DataVector(), maxRank(0), numBins(0), numSeqs(0) {
	try {
        groupmap = NULL; countTable = NULL;
		//set up groupmap for later.
        if (m->groupMode == "group") {
            groupmap = new GroupMap(m->getGroupFile());
            groupmap->readMap(); 
        }else {
            countTable = new CountTable();
            countTable->readTable(m->getCountTableFile(), true, false);
        }

		int hold;
		string inputData;
		f >> label >> hold;
	
		data.assign(hold, "");
		
		for(int i=0;i<hold;i++){
			f >> inputData;
			set(i, inputData);
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "SharedListVector", "SharedListVector");
		exit(1);
	}
}

/***********************************************************************/
void SharedListVector::set(int binNumber, string seqNames){
	try {
		int nNames_old = m->getNumNames(data[binNumber]);
		data[binNumber] = seqNames;
		int nNames_new = m->getNumNames(seqNames);
	
		if(nNames_old == 0)			{	numBins++;				}
		if(nNames_new == 0)			{	numBins--;				}
		if(nNames_new > maxRank)	{	maxRank = nNames_new;	}
	
		numSeqs += (nNames_new - nNames_old);
		
			 
	}
	catch(exception& e) {
		m->errorOut(e, "SharedListVector", "set");
		exit(1);
	}
}

/***********************************************************************/

string SharedListVector::get(int index){
	return data[index];
}

/***********************************************************************/

void SharedListVector::push_back(string seqNames){
	try {
		data.push_back(seqNames);
		int nNames = m->getNumNames(seqNames);
	
		numBins++;
	
		if(nNames > maxRank)	{	maxRank = nNames;	}
	
		numSeqs += nNames;
	}
	catch(exception& e) {
		m->errorOut(e, "SharedListVector", "push_back");
		exit(1);
	}
}

/***********************************************************************/

void SharedListVector::resize(int size){
	data.resize(size);		
}

/***********************************************************************/

int SharedListVector::size(){
	return data.size();
}
/***********************************************************************/

void SharedListVector::clear(){
	numBins = 0;
	maxRank = 0;
	numSeqs = 0;
	return data.clear();
	
}

/***********************************************************************/

void SharedListVector::print(ostream& output){
	try {
		output << label << '\t' << numBins << '\t';
	
		for(int i=0;i<data.size();i++){
			if(data[i] != ""){
				output << data[i] << '\t';
			}
		}
		output << endl;
	}
	catch(exception& e) {
		m->errorOut(e, "SharedListVector", "print");
		exit(1);
	}
}


/***********************************************************************/

RAbundVector SharedListVector::getRAbundVector(){
	try {
		RAbundVector rav;
	
		for(int i=0;i<data.size();i++){
			int binSize = m->getNumNames(data[i]);
			rav.push_back(binSize);
		}
	
	//  This was here before to output data in a nice format, but it screws up the name mapping steps
	//	sort(rav.rbegin(), rav.rend());
	//	
	//	for(int i=data.size()-1;i>=0;i--){
	//		if(rav.get(i) == 0){	rav.pop_back();	}
	//		else{
	//			break;
	//		}
	//	}
		rav.setLabel(label);
	
		return rav;
	}
	catch(exception& e) {
		m->errorOut(e, "SharedListVector", "getRAbundVector");
		exit(1);
	}
}

/***********************************************************************/

SAbundVector SharedListVector::getSAbundVector(){
	try {
		SAbundVector sav(maxRank+1);
	
		for(int i=0;i<data.size();i++){
			int binSize = m->getNumNames(data[i]);	
			sav.set(binSize, sav.get(binSize) + 1);	
		}
		sav.set(0, 0);
		sav.setLabel(label);
	
		return sav;
	}
	catch(exception& e) {
		m->errorOut(e, "SharedListVector", "getSAbundVector");
		exit(1);
	}
}

/***********************************************************************/
SharedOrderVector* SharedListVector::getSharedOrderVector(){
	try {
		SharedOrderVector* order = new SharedOrderVector();
		order->setLabel(label);
	
		for(int i=0;i<numBins;i++){
			int binSize = m->getNumNames(get(i));	//find number of individual in given bin	
			string names = get(i);
            vector<string> binNames;
            m->splitAtComma(names, binNames);
            if (m->groupMode != "group") {
                binSize = 0;
                for (int j = 0; j < binNames.size(); j++) {  binSize += countTable->getNumSeqs(binNames[i]);  }
            }
			for (int j = 0; j < binNames.size(); j++) { 
                if (m->control_pressed) { return order; }
                if (m->groupMode == "group") {
                    string groupName = groupmap->getGroup(binNames[i]);
                    if(groupName == "not found") {	m->mothurOut("Error: Sequence '" + binNames[i] + "' was not found in the group file, please correct."); m->mothurOutEndLine();  exit(1); }
				
                    order->push_back(i, binSize, groupName);  //i represents what bin you are in
                }else {
                    vector<int> groupAbundances = countTable->getGroupCounts(binNames[i]);
                    vector<string> groupNames = countTable->getNamesOfGroups();
                    for (int k = 0; k < groupAbundances.size(); k++) { //groupAbundances.size() == 0 if there is a file mismatch and m->control_pressed is true.
                        if (m->control_pressed) { return order; }
                        for (int l = 0; l < groupAbundances[k]; l++) {  order->push_back(i, binSize, groupNames[k]);  }
                    }
                }
			}
		}

		random_shuffle(order->begin(), order->end());
		order->updateStats();
		
		return order;
	}
	catch(exception& e) {
		m->errorOut(e, "SharedListVector", "getSharedOrderVector");
		exit(1);
	}
}
/***********************************************************************/
SharedRAbundVector SharedListVector::getSharedRAbundVector(string groupName) {
	try {
		SharedRAbundVector rav(data.size());
		
		for(int i=0;i<numBins;i++){
			string names = get(i);
            vector<string> binNames;
            m->splitAtComma(names, binNames);
            for (int j = 0; j < binNames.size(); j++) { 
				if (m->control_pressed) { return rav; }
                if (m->groupMode == "group") {
                    string group = groupmap->getGroup(binNames[j]);
                    if(group == "not found") {	m->mothurOut("Error: Sequence '" + binNames[j] + "' was not found in the group file, please correct."); m->mothurOutEndLine();  exit(1); }
                    if (group == groupName) { //this name is in the group you want the vector for.
                        rav.set(i, rav.getAbundance(i) + 1, group);  //i represents what bin you are in
                    }
                }else {
                    int count = countTable->getGroupCount(binNames[j], groupName);
                    rav.set(i, rav.getAbundance(i) + count, groupName);
                }
			}
		}
		
		rav.setLabel(label);
		rav.setGroup(groupName);

		return rav;
		
	}
	catch(exception& e) {
		m->errorOut(e, "SharedListVector", "getSharedRAbundVector");
		exit(1);
	}
}
/***********************************************************************/
vector<SharedRAbundVector*> SharedListVector::getSharedRAbundVector() {
	try {
		SharedUtil* util;
		util = new SharedUtil();
		vector<SharedRAbundVector*> lookup;  //contains just the groups the user selected
        vector<SharedRAbundVector*> lookupDelete;
		map<string, SharedRAbundVector*> finder;  //contains all groups in groupmap
		
		vector<string> Groups = m->getGroups();
        vector<string> allGroups;
		if (m->groupMode == "group") {  allGroups = groupmap->getNamesOfGroups();  }
        else {  allGroups = countTable->getNamesOfGroups();  }
		util->setGroups(Groups, allGroups);
		m->setGroups(Groups);
		delete util;

		for (int i = 0; i < allGroups.size(); i++) {
			SharedRAbundVector* temp = new SharedRAbundVector(data.size());
			finder[allGroups[i]] = temp;
			finder[allGroups[i]]->setLabel(label);
			finder[allGroups[i]]->setGroup(allGroups[i]);
			if (m->inUsersGroups(allGroups[i], m->getGroups())) {  //if this group is in user groups
				lookup.push_back(finder[allGroups[i]]);
			}else {
                lookupDelete.push_back(finder[allGroups[i]]);
            }
		}
	
		//fill vectors
		for(int i=0;i<numBins;i++){
			string names = get(i);  
			vector<string> binNames;
            m->splitAtComma(names, binNames);
            for (int j = 0; j < binNames.size(); j++) { 
                if (m->groupMode == "group") {
                    string group = groupmap->getGroup(binNames[j]);
                    if(group == "not found") {	m->mothurOut("Error: Sequence '" + binNames[j] + "' was not found in the group file, please correct."); m->mothurOutEndLine();  exit(1); }
                    finder[group]->set(i, finder[group]->getAbundance(i) + 1, group);  //i represents what bin you are in	
                }else{
                    vector<int> counts = countTable->getGroupCounts(binNames[j]);
                    for (int k = 0; k < allGroups.size(); k++) {
                        finder[allGroups[k]]->set(i, finder[allGroups[k]]->getAbundance(i) + counts[k], allGroups[k]);
                    }
                }
			}
		}
        
        for (int j = 0; j < lookupDelete.size(); j++) {  delete lookupDelete[j];  }

		return lookup;
	}
	catch(exception& e) {
		m->errorOut(e, "SharedListVector", "getSharedRAbundVector");
		exit(1);
	}
}

/***********************************************************************/
SharedSAbundVector SharedListVector::getSharedSAbundVector(string groupName) {
	try { 
		SharedSAbundVector sav;
		SharedRAbundVector rav;
		
		rav = this->getSharedRAbundVector(groupName);
		sav = rav.getSharedSAbundVector();
		
		return sav;
	}
	catch(exception& e) {
		m->errorOut(e, "SharedListVector", "getSharedSAbundVector");
		exit(1);
	}
}
/***********************************************************************/

OrderVector SharedListVector::getOrderVector(map<string,int>* orderMap = NULL){
	
	try {
		if(orderMap == NULL){
			OrderVector ov;
		
			for(int i=0;i<data.size();i++){
                string names = data[i];
                vector<string> binNames;
                m->splitAtComma(names, binNames);
				int binSize = binNames.size();	
                if (m->groupMode != "group") {
                    binSize = 0;
                    for (int j = 0; j < binNames.size(); j++) {  binSize += countTable->getNumSeqs(binNames[i]);  }
                }
				for(int j=0;j<binSize;j++){
					ov.push_back(i);
				}
			}
			random_shuffle(ov.begin(), ov.end());
			ov.setLabel(label);
			ov.getNumBins();
		
			return ov;
		
		}
		else{
			OrderVector ov(numSeqs);
		
			for(int i=0;i<data.size();i++){
				string listOTU = data[i];
				vector<string> binNames;
                m->splitAtComma(listOTU, binNames);
                for (int j = 0; j < binNames.size(); j++) { 
                    if(orderMap->count(binNames[j]) == 0){
                        m->mothurOut(binNames[j] + " not found, check *.names file\n");
                        exit(1);
                    }
                    ov.set((*orderMap)[binNames[j]], i);
				}
			}
		
			ov.setLabel(label);
			ov.getNumBins();
		
			return ov;		
		}
	}
	catch(exception& e) {
		m->errorOut(e, "SharedListVector", "getOrderVector");
		exit(1);
	}
}

/***********************************************************************/

