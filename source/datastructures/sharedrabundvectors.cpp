//
//  sharedrabundvectors.cpp
//  Mothur
//
//  Created by Sarah Westcott on 5/15/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#include "sharedrabundvectors.hpp"
#include "sharedutilities.h"

/***********************************************************************/
//reads a shared file
SharedRAbundVectors::SharedRAbundVectors(ifstream& f) : DataVector(){
    try {
        m->clearAllGroups();
        vector<string> allGroups;
        
        int num, count;
        count = 0;
        string holdLabel, nextLabel, groupN;
        
        for (int i = 0; i < lookup.size(); i++) {  if (lookup[i] != NULL) { delete lookup[i];  lookup[i] = NULL; } }  lookup.clear();
        
        //are we at the beginning of the file??
        if (m->saveNextLabel == "") {
            f >> label;
            
            //is this a shared file that has headers
            if (label == "label") {
                //gets "group"
                f >> label; m->gobble(f);
                
                //gets "numOtus"
                f >> label; m->gobble(f);
                
                //eat rest of line
                label = m->getline(f); m->gobble(f);
                
                //parse labels to save
                istringstream iStringStream(label);
                m->sharedBinLabelsInFile.clear();
                while(!iStringStream.eof()){
                    if (m->control_pressed) { break; }
                    string temp;
                    iStringStream >> temp;  m->gobble(iStringStream);
                    
                    m->sharedBinLabelsInFile.push_back(temp);
                }
                
                f >> label >> groupN >> num;
            }else {
                //read in first row since you know there is at least 1 group.
                f >> groupN >> num;
                
                //make binlabels because we don't have any
                string snumBins = toString(num);
                m->sharedBinLabelsInFile.clear();
                for (int i = 0; i < num; i++) {
                    //if there is a bin label use it otherwise make one
                    string binLabel = "Otu";
                    string sbinNumber = toString(i+1);
                    if (sbinNumber.length() < snumBins.length()) {
                        int diff = snumBins.length() - sbinNumber.length();
                        for (int h = 0; h < diff; h++) { binLabel += "0"; }
                    }
                    binLabel += sbinNumber;
                    m->sharedBinLabelsInFile.push_back(binLabel);
                }
            }
        }else {
            label = m->saveNextLabel;
            
            //read in first row since you know there is at least 1 group.
            f >> groupN >> num;
        }
        
        //reset labels, currentLabels may have gotten changed as otus were eliminated because of group choices or sampling
        m->currentSharedBinLabels = m->sharedBinLabelsInFile;
        
        holdLabel = label;
        allGroups.push_back(groupN);
        numBins = num;
        
        //add new vector to lookup
        SharedRAbundVector* temp = new SharedRAbundVector(f, label, groupN, numBins); m->gobble(f);
        push_back(temp);
        
        if (!(f.eof())) { f >> nextLabel; }
        
        //read the rest of the groups info in
        while ((nextLabel == holdLabel) && (f.eof() != true)) {
            f >> groupN >> num;
            SharedRAbundVector* temp = new SharedRAbundVector(f, label, groupN, numBins); m->gobble(f);
            push_back(temp);
            allGroups.push_back(groupN);
            
            if (f.eof() != true) { f >> nextLabel; }
        }
        m->saveNextLabel = nextLabel;
        m->setAllGroups(allGroups);
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundVectors", "SharedRAbundVectors");
        exit(1);
    }
}
/***********************************************************************/
void SharedRAbundVectors::print(ostream& output){
    try {
        sort(lookup.begin(), lookup.end(), compareRAbunds);
        for (int i = 0; i < lookup.size(); i++) {
            if (m->control_pressed) { break; }
            lookup[i]->print(output);
        }
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundVectors", "print");
        exit(1);
    }
}
/***********************************************************************/
void SharedRAbundVectors::printHeaders(ostream& output){
    try {
        string snumBins = toString(numBins);
        output << "label\tGroup\tnumOtus";
        
        if (m->sharedHeaderMode == "tax") {
            for (int i = 0; i < numBins; i++) {
                
                //if there is a bin label use it otherwise make one
                string binLabel = "PhyloType";
                string sbinNumber = toString(i+1);
                if (sbinNumber.length() < snumBins.length()) {
                    int diff = snumBins.length() - sbinNumber.length();
                    for (int h = 0; h < diff; h++) { binLabel += "0"; }
                }
                binLabel += sbinNumber;
                if (i < m->currentSharedBinLabels.size()) {  binLabel = m->currentSharedBinLabels[i]; }
                
                output << '\t' << binLabel ;
            }
            output << endl;
        }else {
            for (int i = 0; i < numBins; i++) {
                //if there is a bin label use it otherwise make one
                string binLabel = "Otu";
                
                if (i < m->currentSharedBinLabels.size()) {  binLabel = m->currentSharedBinLabels[i]; }
                else {
                    string sbinNumber = toString(i+1);
                    if (sbinNumber.length() < snumBins.length()) {
                        int diff = snumBins.length() - sbinNumber.length();
                        for (int h = 0; h < diff; h++) { binLabel += "0"; }
                    }
                    binLabel += sbinNumber;
                }
                
                output  << '\t' << binLabel;
            }
            
            output << endl;
        }
        m->printedSharedHeaders = true;
    }
    catch(exception& e) {
        m->errorOut(e, "SharedVector", "printHeaders");
        exit(1);
    }
}
/***********************************************************************/
int SharedRAbundVectors::push_back(SharedRAbundVector* thisLookup){
    try {
        
        if (numBins == 0) {  numBins = thisLookup->getNumBins();  }
        lookup.push_back(thisLookup);
        sort(lookup.begin(), lookup.end(), compareRAbunds);
        if (label == "") { label = thisLookup->getLabel(); }
        groupNames.clear();
        for (int i = 0; i < lookup.size(); i ++) { groupNames[lookup[i]->getGroup()] = i; }
        
        return lookup.size();
        
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundVectors", "push_back");
        exit(1);
    }
}
/***********************************************************************/
int SharedRAbundVectors::push_back(vector<int> abunds, string binLabel){
    try {
        if (abunds.size() != lookup.size()) {  m->mothurOut("[ERROR]: you have provided " + toString(abunds.size()) + " abundances, but mothur was expecting " + toString(lookup.size()) + ", please correct.\n"); m->control_pressed  = true; return 0; }
        
        for (int i = 0; i < lookup.size(); i ++) { lookup[i]->push_back(abunds[i]); }
        
        if (binLabel == "") { //create one
            map<string, int>::iterator it;
            int otuNum = 0; bool notDone = true;
            
            //find label prefix
            string prefix = "Otu";
            if (m->currentSharedBinLabels[m->currentSharedBinLabels.size()-1][0] == 'P') { prefix = "PhyloType"; }
            
            string tempLabel = m->currentSharedBinLabels[m->currentSharedBinLabels.size()-1];
            string simpleLastLabel = m->getSimpleLabel(tempLabel);
            m->mothurConvert(simpleLastLabel, otuNum); otuNum++;
            string potentialLabel = toString(otuNum);
            
            while (notDone) {
                if (m->control_pressed) { notDone = false; break; }
                
                potentialLabel = toString(otuNum);
                vector<string>::iterator it = find(m->currentSharedBinLabels.begin(), m->currentSharedBinLabels.end(), potentialLabel);
                if (it == m->currentSharedBinLabels.end()) {
                    potentialLabel = prefix + toString(otuNum);
                    it = find(m->currentSharedBinLabels.begin(), m->currentSharedBinLabels.end(), potentialLabel);
                    if (it == m->currentSharedBinLabels.end()) {
                        notDone = false; break;
                    }
                }
                otuNum++;
            }
            
            binLabel = potentialLabel;
        }
        m->sharedBinLabelsInFile.push_back(binLabel);
        m->currentSharedBinLabels.push_back(binLabel);
        
        return lookup.size();
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundVectors", "push_back");
        exit(1);
    }
}
/***********************************************************************/
int SharedRAbundVectors::getOTUTotal(int bin){
    try {
        int totalOTUAbund = 0;
        for (int i = 0; i < lookup.size(); i++) { totalOTUAbund += lookup[i]->get(bin); }
        return totalOTUAbund;
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundVectors", "push_back");
        exit(1);
    }
}
/***********************************************************************/
vector<int> SharedRAbundVectors::getOTU(int bin){
    try {
        vector<int> abunds;
        for (int i = 0; i < lookup.size(); i++) { abunds.push_back(lookup[i]->get(bin)); }
        return abunds;
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundVectors", "push_back");
        exit(1);
    }
}
/***********************************************************************/
void SharedRAbundVectors::setLabels(string l){
    try {
        label = l;
        for (int i = 0; i < lookup.size(); i++) { lookup[i]->setLabel(l); }
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundVectors", "setLabels");
        exit(1);
    }
}
/***********************************************************************/
int SharedRAbundVectors::get(int bin, string group){
    try {
        int abund = 0;
        map<string, int>::iterator it = groupNames.find(group);
        if (it == groupNames.end()) { m->mothurOut("[ERROR]: can not find group " + group + ".\n"); m->control_pressed = true;  }
        else { abund = lookup[it->second]->get(bin); }
        
        return abund;
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundVectors", "get");
        exit(1);
    }
}
/***********************************************************************/
int SharedRAbundVectors::getNumSeqs(string group){
    try {
        int numSeqs = 0;
        map<string, int>::iterator it = groupNames.find(group);
        if (it == groupNames.end()) { m->mothurOut("[ERROR]: can not find group " + group + ".\n"); m->control_pressed = true;  }
        else { numSeqs = lookup[it->second]->getNumSeqs(); }
        
        return numSeqs;
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundVectors", "getNumSeqs");
        exit(1);
    }
}
/***********************************************************************/
void SharedRAbundVectors::set(int bin, int binSize, string group){
    try {
        map<string, int>::iterator it = groupNames.find(group);
        if (it == groupNames.end()) { m->mothurOut("[ERROR]: can not find group " + group + ".\n"); m->control_pressed = true;  }
        else { lookup[it->second]->set(bin, binSize); }
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundVectors", "set");
        exit(1);
    }
}
/***********************************************************************/
int SharedRAbundVectors::removeOTU(int bin){
    try {
        int totalOTUAbund = 0;
        for (int i = 0; i < lookup.size(); i ++) {
            totalOTUAbund += lookup[i]->remove(bin);            
        }
        m->currentSharedBinLabels.erase(m->currentSharedBinLabels.begin()+bin);
        numBins--;
        
        return totalOTUAbund;
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundVectors", "removeOTU");
        exit(1);
    }
}
/***********************************************************************/
vector<string> SharedRAbundVectors::getNamesGroups(){
    try {
        vector<string> names; names.clear();
        for (int i = 0; i < lookup.size(); i ++) { names.push_back(lookup[i]->getGroup()); }
        return names;
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundVectors", "getNamesGroups");
        exit(1);
    }
}
/***********************************************************************/
SharedOrderVector SharedRAbundVectors::getSharedOrderVector(){
    try {
        SharedUtil util;
        
        vector<string> Groups = m->getGroups();
        vector<string> allGroups = m->getAllGroups();
        util.setGroups(Groups, allGroups);
        m->setGroups(Groups);
        
        bool remove = false;
        for (vector<SharedRAbundVector*>::iterator it = lookup.begin(); it != lookup.end();) {
            //if this sharedrabund is not from a group the user wants then delete it.
            if (util.isValidGroup((*it)->getGroup(), m->getGroups()) == false) {
                remove = true;
                delete (*it); (*it) = NULL;
                it = lookup.erase(it);
            }else { ++it; }
        }
        
        if (remove) { eliminateZeroOTUS(); }
        
        SharedOrderVector order;
        for (int i = 0; i < lookup.size(); i++) {
            for (int j = 0; j < lookup[i]->getNumBins(); j++) {
                int abund = lookup[i]->get(j);
                int bin = j+1;
                if (abund != 0) {
                    for (int k = 0; k < abund; k++) {  order.push_back(j, j, lookup[i]->getGroup());  }
                }
            }
        }
        
        return order;
        
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundVector", "getSharedOrderVector");
        exit(1);
    }
}
/***********************************************************************/
void SharedRAbundVectors::removeGroups(vector<string> g){
    try {
        SharedUtil util;
        
        vector<string> Groups = m->getGroups();
        vector<string> allGroups = m->getAllGroups();
        util.setGroups(Groups, allGroups);
        m->setGroups(Groups);
        
        bool remove = false;
        for (vector<SharedRAbundVector*>::iterator it = lookup.begin(); it != lookup.end();) {
            //if this sharedrabund is not from a group the user wants then delete it.
            if (util.isValidGroup((*it)->getGroup(), g) == false) {
                remove = true;
                delete (*it); (*it) = NULL;
                it = lookup.erase(it);
            }else { ++it; }
        }
        
        if (remove) { eliminateZeroOTUS(); }
        
        groupNames.clear();
        for (int i = 0; i < lookup.size(); i ++) { groupNames[lookup[i]->getGroup()] = i; }
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundVector", "removeGroups");
        exit(1);
    }
}
/***********************************************************************/
int SharedRAbundVectors::removeGroups(int minSize, bool silent){
    try {
        vector<string> Groups;
        bool remove = false;
        for (vector<SharedRAbundVector*>::iterator it = lookup.begin(); it != lookup.end();) {
            if ((*it)->getNumSeqs() < minSize) {
                if (!silent) { m->mothurOut((*it)->getGroup() + " contains " + toString((*it)->getNumSeqs()) + ". Eliminating."); m->mothurOutEndLine(); }
                delete (*it); (*it) = NULL;
                it = lookup.erase(it);
            }else {
                Groups.push_back((*it)->getGroup());
                ++it;
            }
        }
        m->setGroups(Groups);
        
        if (remove) { eliminateZeroOTUS(); }
        
        groupNames.clear();
        for (int i = 0; i < lookup.size(); i ++) { groupNames[lookup[i]->getGroup()] = i; }
        
        return lookup.size();
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundVector", "removeGroups");
        exit(1);
    }
}
/***********************************************************************/
int SharedRAbundVectors::getNumSeqsSmallestGroup(){
    try {
        int smallest = 1e6;
        for (int i = 0; i < lookup.size(); i++) {
            if (m->debug) { m->mothurOut("[DEBUG]: " + lookup[i]->getGroup() + " numSeqs = " + toString(lookup[i]->getNumSeqs()) + "\n"); }
            if (lookup[i]->getNumSeqs() < smallest) { smallest = lookup[i]->getNumSeqs(); }
        }
        return smallest;
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundVector", "getNumSeqsSmallestGroup");
        exit(1);
    }
}

/***********************************************************************/
vector<SharedRAbundVector*> SharedRAbundVectors::getSharedRAbundVectors(){
    try {
        SharedUtil util;
        
        vector<string> Groups = m->getGroups();
        vector<string> allGroups = m->getAllGroups();
        util.setGroups(Groups, allGroups);
        m->setGroups(Groups);
        
        vector<SharedRAbundVector*> newLookup;
        bool remove = false;
        for (vector<SharedRAbundVector*>::iterator it = lookup.begin(); it != lookup.end();) {
            if (m->control_pressed) { return newLookup; }
            //if this sharedrabund is not from a group the user wants then delete it.
            if (util.isValidGroup((*it)->getGroup(), m->getGroups()) == false) {
                remove = true;
                delete (*it); (*it) = NULL;
                it = lookup.erase(it);
            }else { ++it; }
        }

        if (remove) { eliminateZeroOTUS(); }
        
        for (int i = 0; i < lookup.size(); i++) {
            if (m->control_pressed) { return newLookup; }
            SharedRAbundVector* temp = new SharedRAbundVector(*lookup[i]);
            newLookup.push_back(temp);
        }
        
        return newLookup;
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundVector", "getSharedRAbundVectors");
        exit(1);
    }
}
/***********************************************************************/
vector<SharedRAbundFloatVector*> SharedRAbundVectors::getSharedRAbundFloatVectors(){
    try {
        SharedUtil util;
        
        vector<string> Groups = m->getGroups();
        vector<string> allGroups = m->getAllGroups();
        util.setGroups(Groups, allGroups);
        m->setGroups(Groups);
        
        vector<SharedRAbundFloatVector*> newLookup;
        bool remove = false;
        for (vector<SharedRAbundVector*>::iterator it = lookup.begin(); it != lookup.end();) {
            if (m->control_pressed) { return newLookup; }
            //if this sharedrabund is not from a group the user wants then delete it.
            if (util.isValidGroup((*it)->getGroup(), m->getGroups()) == false) {
                remove = true;
                delete (*it); (*it) = NULL;
                it = lookup.erase(it);
            }else { ++it; }
        }
        
        if (remove) { eliminateZeroOTUS(); }
        
        for (int i = 0; i < lookup.size(); i++) {
            if (m->control_pressed) { return newLookup; }
            vector<float> abunds;
            vector<int> data = lookup[i]->get();
            string group = lookup[i]->getGroup();
            for (int j = 0; j < data.size(); j++) { abunds.push_back((float)data[j]); }
            SharedRAbundFloatVector* temp = new SharedRAbundFloatVector(abunds);
            temp->setLabel(label);
            temp->setGroup(group);
            newLookup.push_back(temp);
        }
        
        return newLookup;
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundVector", "getSharedRAbundVectors");
        exit(1);
    }
}
/***********************************************************************/
RAbundVector SharedRAbundVectors::getRAbundVector(){
    try {
        SharedUtil util;
        
        vector<string> Groups = m->getGroups();
        vector<string> allGroups = m->getAllGroups();
        util.setGroups(Groups, allGroups);
        m->setGroups(Groups);
        
        bool remove = false;
        for (vector<SharedRAbundVector*>::iterator it = lookup.begin(); it != lookup.end();) {
            //if this sharedrabund is not from a group the user wants then delete it.
            if (util.isValidGroup((*it)->getGroup(), m->getGroups()) == false) {
                remove = true;
                delete (*it); (*it) = NULL;
                it = lookup.erase(it);
            }else { ++it; }
        }
        
        if (remove) { eliminateZeroOTUS(); }
        
        RAbundVector rav;
        for (int i = 0; i < numBins; i++) {
            int abund = getOTUTotal(i);
            rav.push_back(abund);
        }
        
        return rav;
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundVector", "getSharedRAbundVectors");
        exit(1);
    }
}
/***********************************************************************/
RAbundVector SharedRAbundVectors::getRAbundVector(string group){
    try {
        RAbundVector rav;
        rav.setLabel(label);
        
        for (vector<SharedRAbundVector*>::iterator it = lookup.begin(); it != lookup.end();) {
            //if this sharedrabund is not from a group the user wants then delete it.
            if ((*it)->getGroup() == group) {
                for (int i = 0; i < (*it)->getNumBins(); i++) { rav.push_back((*it)->get(i)); }
            }else { ++it; }
        }
        
        return rav;
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundVector", "getRAbundVector");
        exit(1);
    }
}
/***********************************************************************/
SAbundVector SharedRAbundVectors::getSAbundVector(){
    try {
        RAbundVector rav = getRAbundVector();
        return rav.getSAbundVector();
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundVector", "getSharedRAbundVectors");
        exit(1);
    }
}
/***********************************************************************
OrderVector SharedRAbundVectors::getOrderVector(map<string,int>* nameMap = NULL){
    try {
        RAbundVector rav = getRAbundVector();
        return rav.getOrderVector(NULL);
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundVector", "getSharedRAbundVectors");
        exit(1);
    }
}
/**********************************************************************************************************************/
void SharedRAbundVectors::eliminateZeroOTUS() {
    try {
        
        for (int i = 0; i < lookup[0]->getNumBins();) {
            if (m->control_pressed) { break; }
            
            int total = getOTUTotal(i);
            
            //if they are not all zero add this bin
            if (total == 0) { removeOTU(i);  }
            else { i++;  }
        }
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundVectors", "eliminateZeroOTUS");
        exit(1);
    }
}
/***********************************************************************/


