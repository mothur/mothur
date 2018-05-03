//
//  sharedrabundfloatvectors.cpp
//  Mothur
//
//  Created by Sarah Westcott on 5/15/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#include "sharedrabundfloatvectors.hpp"


/***********************************************************************/
//reads a shared file
SharedRAbundFloatVectors::SharedRAbundFloatVectors(ifstream& f, vector<string>& userGroups, string& nextLabel, string& labelTag) : DataVector() {
    try {
        int num;
        string holdLabel, nextLabel, groupN;
        int numUserGroups = userGroups.size();
        printSharedHeaders = true;
        
        for (int i = 0; i < lookup.size(); i++) {  if (lookup[i] != NULL) { delete lookup[i];  lookup[i] = NULL; } }  lookup.clear();
        
        //are we at the beginning of the file??
        if (nextLabel == "") {
            f >> label;
            
            //is this a shared file that has headers
            if (label == "label") {
                //gets "group"
                f >> label; util.gobble(f);
                
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
                    
                    currentLabels.push_back(temp);
                }
                
                if (currentLabels.size() != 0) {
                    string binLabelTag = currentLabels[0];
                    labelTag = "";
                    for (int i = 0; i < binLabelTag.length(); i++) { if (isalpha(binLabelTag[i])){ labelTag += binLabelTag[i]; } }
                }
                f >> label >> groupN >> num;
            }else {
                //read in first row since you know there is at least 1 group.
                f >> groupN >> num;
                
                //make binlabels because we don't have any
                string snumBins = toString(num);
                if (labelTag == "") { labelTag = "Otu"; }
                for (int i = 0; i < num; i++) {
                    //if there is a bin label use it otherwise make one
                    string binLabel = labelTag;
                    string sbinNumber = toString(i+1);
                    if (sbinNumber.length() < snumBins.length()) {
                        int diff = snumBins.length() - sbinNumber.length();
                        for (int h = 0; h < diff; h++) { binLabel += "0"; }
                    }
                    binLabel += sbinNumber;
                    currentLabels.push_back(binLabel);
                }
            }
        }else {
            label = nextLabel;
            
            //read in first row since you know there is at least 1 group.
            f >> groupN >> num;
        }
        
        bool readData = false;
         bool remove = false;
        if (numUserGroups == 0) { //user has not specified groups, so we will use all of them
            userGroups.push_back(groupN);
            readData = true;
        }else{
            if (util.inUsersGroups(groupN, userGroups)) { readData = true; }
            else { remove = true; }// skipline because you are a group we dont care about
        }

        
        holdLabel = label;
        numBins = num;
        
        if (readData) {
            //add new vector to lookup
            SharedRAbundFloatVector* temp = new SharedRAbundFloatVector(f, label, groupN, numBins);
            push_back(temp);
        } else { util.getline(f); }
        util.gobble(f);
        
        if (!(f.eof())) { f >> nextLabel; }
        
        //read the rest of the groups info in
        while ((nextLabel == holdLabel) && (f.eof() != true)) {
            f >> groupN >> num;
            bool readData = false;
            if (numUserGroups == 0) { //user has not specified groups, so we will use all of them
                userGroups.push_back(groupN);
                readData = true;
            }else{
                if (util.inUsersGroups(groupN, userGroups)) { readData = true; }
                else { remove = true; }// skipline because you are a group we dont care about
            }
            
            if (readData) {
                SharedRAbundFloatVector* temp = new SharedRAbundFloatVector(f, label, groupN, numBins);
                push_back(temp);
            }else { util.getline(f); }
            util.gobble(f);
            
            if (f.eof() != true) { f >> nextLabel; }
        }
        if (remove) { eliminateZeroOTUS(); }
        
        otuTag = labelTag;
        
        //error in names of user inputted Groups
        if (lookup.size() < userGroups.size()) { m->mothurOut("[ERROR]: requesting groups not present in files, aborting.\n"); m->setControl_pressed(true); }
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundFloatVectors", "SharedRAbundFloatVectors");
        exit(1);
    }
}
/***********************************************************************/
void SharedRAbundFloatVectors::print(ostream& output){
    try {
        printHeaders(output);
        sort(lookup.begin(), lookup.end(), compareRAbundFloats);
        for (int i = 0; i < lookup.size(); i++) {
            if (m->getControl_pressed()) { break; }
            lookup[i]->print(output);
        }
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundFloatVectors", "print");
        exit(1);
    }
}
/***********************************************************************/
string SharedRAbundFloatVectors::getOTUName(int bin){
    try {
        if (currentLabels.size() < bin) {  }
        else { getOTUNames(); }
        return currentLabels[bin];
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundFloatVectors", "getOTUName");
        exit(1);
    }
}
/***********************************************************************/
void SharedRAbundFloatVectors::setOTUName(int bin, string otuName){
    try {
        if (currentLabels.size() < bin) {  currentLabels[bin] = otuName; }
        else {
            getOTUNames(); //fills currentLabels if needed
            if (currentLabels.size() < bin) {  currentLabels[bin] = otuName; }
            else {
                m->setControl_pressed(true);
                m->mothurOut("[ERROR]: " + toString(bin) + " bin does not exist\n");
            }
        }
        
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundFloatVectors", "setOTUName");
        exit(1);
    }
}
/***********************************************************************/
int SharedRAbundFloatVectors::push_back(SharedRAbundFloatVector* thisLookup){
    try {
        if (numBins == 0) { numBins = thisLookup->getNumBins();  }
        lookup.push_back(thisLookup);
        sort(lookup.begin(), lookup.end(), compareRAbundFloats);
        if (label == "") { label = thisLookup->getLabel(); }
        groupNames.clear();
        for (int i = 0; i < lookup.size(); i ++) { groupNames[lookup[i]->getGroup()] = i; }
        return lookup.size();
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundFloatVectors", "push_back");
        exit(1);
    }
}
/***********************************************************************/
float SharedRAbundFloatVectors::getOTUTotal(int bin){
    try {
        float totalOTUAbund = 0;
        for (int i = 0; i < lookup.size(); i++) { totalOTUAbund += lookup[i]->get(bin); }
        return totalOTUAbund;
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundFloatVectors", "getOTUTotal");
        exit(1);
    }
}
/***********************************************************************/
vector<float> SharedRAbundFloatVectors::getOTU(int bin){
    try {
        vector<float> abunds;
        for (int i = 0; i < lookup.size(); i++) { abunds.push_back(lookup[i]->get(bin)); }
        return abunds;
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundFloatVectors", "getOTU");
        exit(1);
    }
}

/***********************************************************************/
void SharedRAbundFloatVectors::setLabels(string l){
    try {
        label = l;
        for (int i = 0; i < lookup.size(); i++) { lookup[i]->setLabel(l); }
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundFloatVectors", "setLabels");
        exit(1);
    }
}
/***********************************************************************/
float SharedRAbundFloatVectors::get(int bin, string group){
    try {
        float abund = 0;
        map<string, int>::iterator it = groupNames.find(group);
        if (it == groupNames.end()) { m->mothurOut("[ERROR]: can not find group " + group + ".\n"); m->setControl_pressed(true);  }
        else { abund = lookup[it->second]->get(bin); }
        return abund;
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundFloatVectors", "get");
        exit(1);
    }
}
/***********************************************************************/
float SharedRAbundFloatVectors::getNumSeqs(string group){
    try {
        float numSeqs = 0;
        map<string, int>::iterator it = groupNames.find(group);
        if (it == groupNames.end()) { m->mothurOut("[ERROR]: can not find group " + group + ".\n"); m->setControl_pressed(true);  }
        else { numSeqs = lookup[it->second]->getNumSeqs(); }
        
        return numSeqs;
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundFloatVectors", "getNumSeqs");
        exit(1);
    }
}
/***********************************************************************/
void SharedRAbundFloatVectors::set(int bin, float binSize, string group){
    try {
        map<string, int>::iterator it = groupNames.find(group);
        if (it == groupNames.end()) { m->mothurOut("[ERROR]: can not find group " + group + ".\n"); m->setControl_pressed(true);  }
        else { lookup[it->second]->set(bin, binSize); }
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundFloatVectors", "set");
        exit(1);
    }
}
/***********************************************************************/
float SharedRAbundFloatVectors::removeOTU(int bin){
    try {
        float totalOTUAbund = 0;
        for (int i = 0; i < lookup.size(); i ++) { totalOTUAbund += lookup[i]->remove(bin); }
        
        currentLabels.erase(currentLabels.begin()+bin);
        numBins--;
        
        return totalOTUAbund;
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundFloatVectors", "removeOTU");
        exit(1);
    }
}
/***********************************************************************/
void SharedRAbundFloatVectors::setOTUNames(vector<string> names){
    try {
        currentLabels.clear();
        currentLabels = names;
        getOTUNames();
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundFloatVectors", "setOTUNames");
        exit(1);
    }
}
/***********************************************************************/
vector<string> SharedRAbundFloatVectors::getOTUNames(){
    try {
        util.getOTUNames(currentLabels, numBins, otuTag);
        return currentLabels;
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundFloatVectors", "getOTUNames");
        exit(1);
    }
}
/***********************************************************************/
void SharedRAbundFloatVectors::printHeaders(ostream& output){
    try {
        if (printSharedHeaders) {
            getOTUNames();
        
            output << "label\tGroup\tnumOtus";
            for (int i = 0; i < numBins; i++) { output  << '\t' << currentLabels[i]; } output << endl;
        
            printSharedHeaders = false;
        }
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundFloatVectors", "printHeaders");
        exit(1);
    }
}
/***********************************************************************/
vector<string> SharedRAbundFloatVectors::getNamesGroups(){
    try {
        vector<string> names;
        for (int i = 0; i < lookup.size(); i ++) { names.push_back(lookup[i]->getGroup()); }
        return names;
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundFloatVectors", "getNamesGroups");
        exit(1);
    }
}
/***********************************************************************/
float SharedRAbundFloatVectors::getNumSeqsSmallestGroup(){
    try {
        float smallest = 1e6;
        for (int i = 0; i < lookup.size(); i++) {
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
RAbundVector SharedRAbundFloatVectors::getRAbundVector(){
    try {
        RAbundVector rav;
        for (int i = 0; i < numBins; i++) {
            float abund = getOTUTotal(i);
            rav.push_back((int)abund);
        }
        
        return rav;
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundFloatVectors", "getSharedRAbundVectors");
        exit(1);
    }
}
/***********************************************************************/
SAbundVector SharedRAbundFloatVectors::getSAbundVector(){
    try {
        RAbundVector rav = getRAbundVector();
        return rav.getSAbundVector();
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundFloatVectors", "getSharedRAbundVectors");
        exit(1);
    }
}
/***********************************************************************/
vector<SharedRAbundFloatVector*> SharedRAbundFloatVectors::getSharedRAbundFloatVectors(){
    try {
        vector<SharedRAbundFloatVector*> newLookup;
        for (int i = 0; i < lookup.size(); i++) {
            SharedRAbundFloatVector* temp = new SharedRAbundFloatVector(*lookup[i]);
            newLookup.push_back(temp);
        }
        
        return newLookup;
        
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundFloatVectors", "getSharedRAbundVectors");
        exit(1);
    }
}
/***********************************************************************/
vector<SharedRAbundVector*> SharedRAbundFloatVectors::getSharedRAbundVectors(){
    try {
        vector<SharedRAbundVector*> newLookup;
        for (int i = 0; i < lookup.size(); i++) {
            SharedRAbundVector* temp = new SharedRAbundVector(lookup[i]->getSharedRAbundVector());
            newLookup.push_back(temp);
        }

        return newLookup;
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundFloatVectors", "getSharedRAbundVectors");
        exit(1);
    }
}
/***********************************************************************/
void SharedRAbundFloatVectors::removeGroups(vector<string> g){
    try {
        bool remove = false;
        for (vector<SharedRAbundFloatVector*>::iterator it = lookup.begin(); it != lookup.end();) {
            //if this sharedrabund is not from a group the user wants then delete it.
            if (util.inUsersGroups((*it)->getGroup(), g) ) {
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
        m->errorOut(e, "SharedRAbundFloatVectors", "removeGroups");
        exit(1);
    }
}
/***********************************************************************/
int SharedRAbundFloatVectors::removeGroups(int minSize, bool silent){
    try {
        vector<string> Groups;
        bool remove = false;
        for (vector<SharedRAbundFloatVector*>::iterator it = lookup.begin(); it != lookup.end();) {
            if ((*it)->getNumSeqs() < minSize) {
                if (!silent) { m->mothurOut((*it)->getGroup() + " contains " + toString((*it)->getNumSeqs()) + ". Eliminating."); m->mothurOutEndLine(); }
                delete (*it); (*it) = NULL;
                it = lookup.erase(it);
            }else {
                Groups.push_back((*it)->getGroup());
                ++it;
            }
        }
        
        if (remove) { eliminateZeroOTUS(); }
        
        groupNames.clear();
        for (int i = 0; i < lookup.size(); i ++) { groupNames[lookup[i]->getGroup()] = i; }
        
        return lookup.size();
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundFloatVectors", "removeGroups");
        exit(1);
    }
}
/**********************************************************************************************************************/
void SharedRAbundFloatVectors::eliminateZeroOTUS() {
    try {
        if (lookup.size() > 1) {
            for (int i = 0; i < lookup[0]->getNumBins();) {
                if (m->getControl_pressed()) { break; }
                
                float total = getOTUTotal(i);
                
                //if they are not all zero add this bin
                if (total == 0) { removeOTU(i); }
                else { ++i;  }
            }
        }
    }
    catch(exception& e) {
        m->errorOut(e, "SharedRAbundFloatVectors", "eliminateZeroOTUS");
        exit(1);
    }
}
/***********************************************************************/

