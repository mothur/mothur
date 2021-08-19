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


/***********************************************************************/
SharedListVector::SharedListVector(ifstream& f, vector<string>& userGroups, string& previousLabel, string& labelTag) : DataVector(), maxRank(0), numBins(0), numSeqs(0) {
	try {
        Utils util;
        groups = userGroups; fillGroups = true;
        if (groups.size() > 0) { fillGroups = false; }
        
        CurrentFile* current = CurrentFile::getInstance();
        groupMode = current->getGroupMode();
        groupmap = NULL; countTable = NULL;
		//set up groupmap for later.
        if (groupMode == "group") {
            groupmap = new GroupMap(current->getGroupFile());
            groupmap->readMap();
            if (fillGroups) { groups = groupmap->getNamesOfGroups(); m->mothurOut("[ERROR]: requesting groups not present in files, aborting.\n"); fillGroups = false; }
            else { if (!util.isSubset(groupmap->getNamesOfGroups(), groups)) { m->mothurOut("[ERROR]: requesting groups not present in files, aborting.\n"); m->setControl_pressed(true); } }
        }else {
            countTable = new CountTable();
            countTable->readTable(current->getCountFile(), true, false);
            if (fillGroups) { groups = countTable->getNamesOfGroups(); fillGroups = false; }
            else { if (!util.isSubset(countTable->getNamesOfGroups(), groups)) { m->mothurOut("[ERROR]: requesting groups not present in files, aborting.\n"); m->setControl_pressed(true); } }
        }

        int hold;
        
		//are we at the beginning of the file?? If yes, read or create headers
		if (previousLabel == "") {
			f >> label;
            
			//is this a shared file that has headers
			if (label == "label") {
				
				//gets "numOtus"
				f >> label; util.gobble(f);
				
				//eat rest of line
				label = util.getline(f); util.gobble(f);
				
				//parse labels to save
				istringstream iStringStream(label);
				while(!iStringStream.eof()){
					if (m->getControl_pressed()) { break; }
					string temp;
					iStringStream >> temp;  util.gobble(iStringStream);
                    
					binLabels.push_back(temp);
				}
                if (binLabels.size() != 0) {
                    string binLabelTag = binLabels[0];
                    labelTag = "";
                    for (int i = 0; i < binLabelTag.length(); i++) { if (isalpha(binLabelTag[i])){ labelTag += binLabelTag[i]; } }
                }
				f >> label >> hold;
			}else {
                //read in first row
                f >> hold;
                
                //make binlabels because we don't have any
                string snumBins = toString(hold);
                if (labelTag == "") { labelTag = "Otu"; }
                for (int i = 0; i < hold; i++) {
                    //if there is a bin label use it otherwise make one
                    string binLabel = labelTag;
                    string sbinNumber = toString(i+1);
                    if (sbinNumber.length() < snumBins.length()) {
                        int diff = snumBins.length() - sbinNumber.length();
                        for (int h = 0; h < diff; h++) { binLabel += "0"; }
                    }
                    binLabel += sbinNumber;
                    binLabels.push_back(binLabel);
                }
            }
            
		}else {  f >> label >> hold; }
        
		data.assign(hold, "");
		string inputData = "";
        otuTag = labelTag;
        previousLabel = label;
        
		for(int i=0;i<hold;i++){
			f >> inputData;
			set(i, inputData);
		}
		util.gobble(f);
	}
	catch(exception& e) {
		m->errorOut(e, "SharedListVector", "SharedListVector");
		exit(1);
	}
}

/***********************************************************************/
void SharedListVector::set(int binNumber, string seqNames){
	try {
        Utils util;
		int nNames_old = util.getNumNames(data[binNumber]);
		data[binNumber] = seqNames;
		int nNames_new = util.getNumNames(seqNames);
	
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

void SharedListVector::setLabels(vector<string> labels){
	try {
		binLabels = labels;
        getLabels();
	}
	catch(exception& e) {
		m->errorOut(e, "SharedListVector", "setLabels");
		exit(1);
	}
}

/***********************************************************************/
//could potentially end up with duplicate binlabel names with code below.
//we don't currently use them in a way that would do that.
//if you had a listfile that had been subsampled and then added to it, dup names would be possible.
vector<string> SharedListVector::getLabels(){
    try {
        Utils util; util.getOTUNames(binLabels, numBins, otuTag);
        return binLabels;
    }
	catch(exception& e) {
		m->errorOut(e, "SharedListVector", "getLabels");
		exit(1);
	}
}
/***********************************************************************/

void SharedListVector::push_back(string seqNames){
	try {
        Utils util;
		data.push_back(seqNames);
		int nNames = util.getNumNames(seqNames);
	
		numBins++;
	
		if(nNames > maxRank)	{	maxRank = nNames;	}
	
		numSeqs += nNames;
    
        int otuNum = numBins; bool notDone = true;
        
        //find label prefix
        string prefix = "Otu";
        if (binLabels[binLabels.size()-1][0] == 'P') { prefix = "PhyloType"; }
        
        string tempLabel = binLabels[binLabels.size()-1];
        string simpleLastLabel = util.getSimpleLabel(tempLabel);
        util.mothurConvert(simpleLastLabel, otuNum); otuNum++;
        string potentialLabel = toString(otuNum);
        
        while (notDone) {
            if (m->getControl_pressed()) { notDone = false; break; }
            
            potentialLabel = toString(otuNum);
            vector<string>::iterator it = find(binLabels.begin(), binLabels.end(), potentialLabel);
            if (it == binLabels.end()) {
                potentialLabel = prefix + toString(otuNum);
                it = find(binLabels.begin(), binLabels.end(), potentialLabel);
                if (it == binLabels.end()) {
                    notDone = false; break;
                }
            }
            otuNum++;
        }
        
       
        binLabels.push_back(potentialLabel);

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
    binLabels.clear();
	return data.clear();
	
}

/***********************************************************************/

void SharedListVector::print(ostream& output){
	try {
		output << label << '\t' << numBins;
	
		for(int i=0;i<data.size();i++){
			if(data[i] != ""){
				output << '\t' << data[i];
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
        Utils util;
		for(int i=0;i<data.size();i++){
			int binSize = util.getNumNames(data[i]);
			rav.push_back(binSize);
		}
	
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
        Utils util;
		for(int i=0;i<data.size();i++){
			int binSize = util.getNumNames(data[i]);	
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
        Utils util;
		for(int i=0;i<numBins;i++){

			string names = get(i);
            vector<string> binNames;
            util.splitAtComma(names, binNames);
            
			for (int j = 0; j < binNames.size(); j++) { 
                if (m->getControl_pressed()) { return order; }
                if (groupMode == "group") {
                    string groupName = groupmap->getGroup(binNames[i]);
                    if(groupName == "not found") {	m->mothurOut("Error: Sequence '" + binNames[i] + "' was not found in the group file, please correct.\n");  exit(1); }
                    if (util.inUsersGroups(groupName, groups)) { order->push_back(i,  groupName);  }//i represents what bin you are in
                }else {
                    vector<int> groupAbundances = countTable->getGroupCounts(binNames[i]);
                    vector<string> groupNames = countTable->getNamesOfGroups();
                    for (int k = 0; k < groupAbundances.size(); k++) { //groupAbundances.size() == 0 if there is a file mismatch and m->control_pressed is true.
                        if (m->getControl_pressed()) { return order; }
                        for (int l = 0; l < groupAbundances[k]; l++) { //for each abundance != 0, add a individual for each
                            if (util.inUsersGroups(groupNames[k], groups)) { order->push_back(i, groupNames[k]); }
                        }
                    }
                }
			}
		}

		util.mothurRandomShuffle(*order);
		order->updateStats();
		
		return order;
	}
	catch(exception& e) {
		m->errorOut(e, "SharedListVector", "getSharedOrderVector");
		exit(1);
	}
}
/***********************************************************************/
SharedRAbundVectors* SharedListVector::getSharedRAbundVector() {
	try {
        vector<SharedRAbundVector*> lookup;  //contains just the groups the user selected
        //vector<string> groups;
        map<string, SharedRAbundVector*> finder;  //contains all groups in groupmap
        map<string, SharedRAbundVector*>::iterator it;
        
		for (int i = 0; i < groups.size(); i++) {
            SharedRAbundVector* temp = new SharedRAbundVector(numBins);
            finder[groups[i]] = temp;
            finder[groups[i]]->setLabel(label);
            finder[groups[i]]->setGroup(groups[i]);
            lookup.push_back(finder[groups[i]]);
        }
        Utils util;
		//fill vectors
		for(int i=0;i<numBins;i++){
			string names = get(i);  
			vector<string> binNames;
            util.splitAtComma(names, binNames);
            for (int j = 0; j < binNames.size(); j++) { 
                if (groupMode == "group") {
                    string group = groupmap->getGroup(binNames[j]);
                    if(group == "not found") {	m->mothurOut("Error: Sequence '" + binNames[j] + "' was not found in the group file, please correct.\n");   exit(1); }
                    it = finder.find(group);
                    if (it != finder.end()) { it->second->set(i, it->second->get(i) + 1); } //i represents what bin you are in
                }else{
                    vector<int> counts = countTable->getGroupCounts(binNames[j]);
                    vector<string> allGroups = countTable->getNamesOfGroups();
                    for (int k = 0; k < allGroups.size(); k++) {
                        it = finder.find(allGroups[k]);
                        if (it != finder.end()) { it->second->set(i, it->second->get(i) + counts[k]); } //i represents what bin you are in
                    }
                }
			}
		}
        
        SharedRAbundVectors* shared = new SharedRAbundVectors(otuTag);
        for (int j = 0; j < lookup.size(); j++) {  shared->push_back(lookup[j]);  }
        shared->setOTUNames(binLabels);
        shared->eliminateZeroOTUS();

		return shared;
	}
	catch(exception& e) {
		m->errorOut(e, "SharedListVector", "getSharedRAbundVector");
		exit(1);
	}
}
/***********************************************************************/
SharedRAbundFloatVectors* SharedListVector::getSharedRAbundFloatVector() {
    try {
        SharedRAbundVectors* shared = getSharedRAbundVector();
        vector<SharedRAbundFloatVector*> thisLookup = shared->getSharedRAbundFloatVectors();
        SharedRAbundFloatVectors* sharedFloat = new SharedRAbundFloatVectors(otuTag);
        for (int j = 0; j < thisLookup.size(); j++) {  sharedFloat->push_back(thisLookup[j]);  }
        return sharedFloat;
    }
    catch(exception& e) {
        m->errorOut(e, "SharedListVector", "getSharedRAbundVector");
        exit(1);
    }
}

/***********************************************************************/

OrderVector SharedListVector::getOrderVector(map<string,int>* orderMap = NULL){
	
	try {
        Utils util;
		if(orderMap == NULL){
			OrderVector ov;
		
			for(int i=0;i<data.size();i++){
                string names = data[i];
                vector<string> binNames;
                util.splitAtComma(names, binNames);
				int binSize = binNames.size();	
                if (groupMode != "group") {
                    binSize = 0;
                    for (int j = 0; j < binNames.size(); j++) {  binSize += countTable->getNumSeqs(binNames[i]);  }
                }
				for(int j=0;j<binSize;j++){
					ov.push_back(i);
				}
			}
			util.mothurRandomShuffle(ov);
			ov.setLabel(label);
			ov.getNumBins();
		
			return ov;
		
		}
		else{
			OrderVector ov(numSeqs);
		
			for(int i=0;i<data.size();i++){
				string listOTU = data[i];
				vector<string> binNames;
                util.splitAtComma(listOTU, binNames);
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

