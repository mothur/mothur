/*
 *  sharedrabundfloatvector.cpp
 *  Mothur
 *
 *  Created by westcott on 8/18/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "sharedrabundfloatvector.h"
#include "sharedutilities.h"

/***********************************************************************/

SharedRAbundFloatVector::SharedRAbundFloatVector() : DataVector(), maxRank(0.0), numBins(0), numSeqs(0.0) {}
/***********************************************************************/

SharedRAbundFloatVector::~SharedRAbundFloatVector() {}

/***********************************************************************/
SharedRAbundFloatVector::SharedRAbundFloatVector(int n) : DataVector(), maxRank(0.0), numBins(n), numSeqs(0.0) {
		
		individualFloat newGuy;
		//initialize data
		for (int i=0; i< n; i++) {
			newGuy.bin = i;
			newGuy.abundance = 0.0;
			data.push_back(newGuy);
		}
}
/***********************************************************************/
//reads a shared file
SharedRAbundFloatVector::SharedRAbundFloatVector(ifstream& f) : DataVector(), maxRank(0.0), numBins(0), numSeqs(0.0) {
	try {
		
		m->clearAllGroups();
		vector<string> allGroups;
		
		int num, count;
		float inputData;
		count = 0;  
		string holdLabel, nextLabel, groupN;
		individualFloat newguy;
		
		for (int i = 0; i < lookup.size(); i++) {  delete lookup[i]; lookup[i] = NULL; }
		lookup.clear();
		
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
		
		//add new vector to lookup
		SharedRAbundFloatVector* temp = new SharedRAbundFloatVector();
		lookup.push_back(temp);
		lookup[0]->setLabel(label);
		lookup[0]->setGroup(groupN);
		
		allGroups.push_back(groupN);
		
		//fill vector.  data = first sharedrabund in file
		for(int i=0;i<num;i++){
			f >> inputData;
			
			lookup[0]->push_back(inputData, groupN); //abundance, bin, group
			push_back(inputData, groupN);
			
			if (inputData > maxRank) { maxRank = inputData; }
		}
		
		m->gobble(f);
		
		if (f.eof() != true) { f >> nextLabel; }
		
		//read the rest of the groups info in
		while ((nextLabel == holdLabel) && (f.eof() != true)) {
			f >> groupN >> num;
            
			count++;
			
			allGroups.push_back(groupN);
			
			//add new vector to lookup
			temp = new SharedRAbundFloatVector();
			lookup.push_back(temp);
			lookup[count]->setLabel(label);
			lookup[count]->setGroup(groupN);

			//fill vector.  
			for(int i=0;i<num;i++){
				f >> inputData;
				lookup[count]->push_back(inputData, groupN); //abundance, bin, group
			}
			
			m->gobble(f);
				
			if (f.eof() != true) { f >> nextLabel; }
		}
		
		m->saveNextLabel = nextLabel;
		m->setAllGroups(allGroups);
	}
	catch(exception& e) {
		m->errorOut(e, "SharedRAbundFloatVector", "SharedRAbundFloatVector");
		exit(1);
	}
}

/***********************************************************************/

void SharedRAbundFloatVector::set(int binNumber, float newBinSize, string groupname){
	try {
		float oldBinSize = data[binNumber].abundance;
		data[binNumber].abundance = newBinSize;
		data[binNumber].group = groupname;
		
		if(newBinSize > maxRank)	{	newBinSize = newBinSize;	}
	
		numSeqs += (newBinSize - oldBinSize);
	}
	catch(exception& e) {
		m->errorOut(e, "SharedRAbundFloatVector", "set");
		exit(1);
	}
}
/***********************************************************************/

void SharedRAbundFloatVector::clear(){
	numBins = 0;
	maxRank = 0;
	numSeqs = 0;
	data.clear();
	for (int i = 0; i < lookup.size(); i++) {  delete lookup[i]; lookup[i] = NULL; }
	lookup.clear();
}
/***********************************************************************/
float SharedRAbundFloatVector::getAbundance(int index){
	return data[index].abundance;	
}
/***********************************************************************/
//returns vector of abundances
vector<float> SharedRAbundFloatVector::getAbundances(){
    vector<float> abunds;
    for (int i = 0; i < data.size(); i++) {
        abunds.push_back(data[i].abundance);
    }
    
	return abunds;
}
/***********************************************************************/
individualFloat SharedRAbundFloatVector::get(int index){
	return data[index];	
}
/***********************************************************************/
void SharedRAbundFloatVector::push_back(float binSize, string groupName){
	try {
		individualFloat newGuy;
		newGuy.abundance = binSize;
		newGuy.group = groupName;
		newGuy.bin = data.size();
		
		data.push_back(newGuy);
		numBins++;
	
		if(binSize > maxRank){	maxRank = binSize;	}
	
		numSeqs += binSize;
	}
	catch(exception& e) {
		m->errorOut(e, "SharedRAbundFloatVector", "push_back");
		exit(1);
	}
}
/***********************************************************************/
void SharedRAbundFloatVector::insert(float binSize, int otu, string groupName){
	try {
		individualFloat newGuy;
		newGuy.abundance = binSize;
		newGuy.group = groupName;
		newGuy.bin = otu;
		
		data.insert(data.begin()+otu, newGuy);
		numBins++;
	
		if(binSize > maxRank){	maxRank = binSize;	}

		numSeqs += binSize;
	}
	catch(exception& e) {
		m->errorOut(e, "SharedRAbundFloatVector", "insert");
		exit(1);
	}
}

/***********************************************************************/
void SharedRAbundFloatVector::push_front(float binSize, int otu, string groupName){
	try {
		individualFloat newGuy;
		newGuy.abundance = binSize;
		newGuy.group = groupName;
		newGuy.bin = otu;
		
		data.insert(data.begin(), newGuy);
		numBins++;
	
		if(binSize > maxRank){	maxRank = binSize;	}
	
		numSeqs += binSize;
	}
	catch(exception& e) {
		m->errorOut(e, "SharedRAbundFloatVector", "push_front");
		exit(1);
	}
}
/**********************************************************************/
void SharedRAbundFloatVector::pop_back(){
	numSeqs -= data[data.size()-1].abundance;
	numBins--;
	data.pop_back();
}
/***********************************************************************/
void SharedRAbundFloatVector::resize(int size){
	data.resize(size);
}
/**********************************************************************/
int SharedRAbundFloatVector::size(){
	return data.size();
}
/***********************************************************************/
void SharedRAbundFloatVector::printHeaders(ostream& output){
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
				
				output  << '\t' << binLabel;
			}
			output << endl;
		}else {
			for (int i = 0; i < numBins; i++) {  
				//if there is a bin label use it otherwise make one
				string binLabel = "Otu";
				string sbinNumber = toString(i+1);
				if (sbinNumber.length() < snumBins.length()) { 
					int diff = snumBins.length() - sbinNumber.length();
					for (int h = 0; h < diff; h++) { binLabel += "0"; }
				}
				binLabel += sbinNumber;
				if (i < m->currentSharedBinLabels.size()) {  binLabel = m->currentSharedBinLabels[i]; }
				
				output  << '\t' << binLabel;
			}
			
			output << endl;
		}
		
		m->printedSharedHeaders = true;
	}
	catch(exception& e) {
		m->errorOut(e, "SharedRAbundVector", "printHeaders");
		exit(1);
	}
}
/***********************************************************************/
void SharedRAbundFloatVector::print(ostream& output){
	try {
		output << numBins;
	
		for(int i=0;i<data.size();i++){		output  << '\t' << data[i].abundance;		}
		output << endl;
	}
	catch(exception& e) {
		m->errorOut(e, "SharedRAbundFloatVector", "print");
		exit(1);
	}
}
/***********************************************************************/
string SharedRAbundFloatVector::getGroup(){
	return group;
}
/***********************************************************************/
void SharedRAbundFloatVector::setGroup(string groupName){
	group = groupName;
}
/***********************************************************************/
int SharedRAbundFloatVector::getGroupIndex()  { return index; }
/***********************************************************************/
void SharedRAbundFloatVector::setGroupIndex(int vIndex)	{ index = vIndex; }
/***********************************************************************/
int SharedRAbundFloatVector::getNumBins(){	return numBins;	}
/***********************************************************************/
float SharedRAbundFloatVector::getNumSeqs(){	return numSeqs;	}
/***********************************************************************/
float SharedRAbundFloatVector::getMaxRank(){	return maxRank;	}
/***********************************************************************/
SharedRAbundFloatVector SharedRAbundFloatVector::getSharedRAbundFloatVector(){
	return *this;			
}
/***********************************************************************
SharedRAbundVector SharedRAbundFloatVector::getSharedRAbundVector(){
	try {
		SharedRAbundVector rav(numBins);
		rav.setLabel(label);
		rav.setGroup(group);
		
		for (int i = 0; i < data.size(); i++) {
			
			rav.push_back(data[i].abundance);
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "SharedRAbundFloatVector", "getSharedRAbundVector");
		exit(1);
	}		
}
***********************************************************************/
vector<SharedRAbundFloatVector*> SharedRAbundFloatVector::getSharedRAbundFloatVectors(){
	try {
		SharedUtil* util;
		util = new SharedUtil();
		
		vector<string> Groups = m->getGroups();
		vector<string> allGroups = m->getAllGroups();
		util->setGroups(Groups, allGroups);
		m->setGroups(Groups);
		
		bool remove = false;
		for (int i = 0; i < lookup.size(); i++) {
			//if this sharedrabund is not from a group the user wants then delete it.
			if (util->isValidGroup(lookup[i]->getGroup(), m->getGroups()) == false) { 
				delete lookup[i]; lookup[i] = NULL;
				lookup.erase(lookup.begin()+i); 
				i--; 
				remove = true;
			}
		}
		
		delete util;
		
		if (remove) { eliminateZeroOTUS(); }
	
		return lookup;
	}
	catch(exception& e) {
		m->errorOut(e, "SharedRAbundFloatVector", "getSharedRAbundFloatVectors");
		exit(1);
	}
}
/***********************************************************************/

RAbundVector SharedRAbundFloatVector::getRAbundVector() {
	try {
		RAbundVector rav(numBins);
		
		//this is not functional, not sure how to handle it yet, but I need the stub because it is a pure function
		
		rav.setLabel(label);
		return rav;
	}
	catch(exception& e) {
		m->errorOut(e, "SharedRAbundFloatVector", "getRAbundVector");
		exit(1);
	}
}
/***********************************************************************

SharedSAbundVector SharedRAbundVector::getSharedSAbundVector(){
	try {
		SharedSAbundVector sav(maxRank+1);
		
		for(int i=0;i<data.size();i++){
			int abund = data[i].abundance;
			sav.set(abund, sav.getAbundance(abund) + 1, group);
		}
		
		sav.set(0, 0, group);
		sav.setLabel(label);
		sav.setGroup(group);
		
		return sav;
	}
	catch(exception& e) {
		m->errorOut(e, "SharedRAbundVector", "getSharedSAbundVector");
		exit(1);
	}
}
***********************************************************************/

SAbundVector SharedRAbundFloatVector::getSAbundVector() {
	try {
		SAbundVector sav(ceil(maxRank)+1);
		
		//this is not functional, not sure how to handle it yet, but I need the stub because it is a pure function
		
		sav.set(0, 0);
		sav.setLabel(label);
		return sav;
	}
	catch(exception& e) {
		m->errorOut(e, "SharedRAbundFloatVector", "getSAbundVector");		
		exit(1);
	}
}

/***********************************************************************

SharedOrderVector SharedRAbundFloatVector::getSharedOrderVector() {
	try {
		SharedOrderVector ov;
	
		for(int i=0;i<data.size();i++){
			int round = ceil(data[i].abundance);
			for(int j=0;j<round;j++){
				ov.push_back(data[i].bin, round, data[i].group);
			}
		}
		random_shuffle(ov.begin(), ov.end());

		ov.setLabel(label);	
		ov.updateStats();
		
		return ov;
	}
	catch(exception& e) {
		m->errorOut(e, "SharedRAbundFloatVector", "getSharedOrderVector");
		exit(1);
	}
}
***********************************************************************/
//this is not functional, not sure how to handle it yet, but I need the stub because it is a pure function
OrderVector SharedRAbundFloatVector::getOrderVector(map<string,int>* nameMap = NULL) {
	try {
		OrderVector ov;
	
		for(int i=0;i<data.size();i++){
			int round = ceil(data[i].abundance);
			for(int j=0;j<round;j++){
				ov.push_back(i);
			}
		}
		random_shuffle(ov.begin(), ov.end());

		ov.setLabel(label);	
		return ov;
	}
	catch(exception& e) {
		m->errorOut(e, "SharedRAbundFloatVector", "getOrderVector");
		exit(1);
	}
}
//**********************************************************************************************************************
int SharedRAbundFloatVector::eliminateZeroOTUS() {
	try {
		
		vector<SharedRAbundFloatVector*> newLookup;
		for (int i = 0; i < lookup.size(); i++) {
			SharedRAbundFloatVector* temp = new SharedRAbundFloatVector();
			temp->setLabel(lookup[i]->getLabel());
			temp->setGroup(lookup[i]->getGroup());
			newLookup.push_back(temp);
		}
		
		//for each bin
		vector<string> newBinLabels;
		string snumBins = toString(lookup[0]->getNumBins());
		for (int i = 0; i < lookup[0]->getNumBins(); i++) {
			if (m->control_pressed) { for (int j = 0; j < newLookup.size(); j++) {  delete newLookup[j];  } return 0; }
			
			//look at each sharedRabund and make sure they are not all zero
			bool allZero = true;
			for (int j = 0; j < lookup.size(); j++) {
				if (lookup[j]->getAbundance(i) != 0) { allZero = false;  break;  }
			}
			
			//if they are not all zero add this bin
			if (!allZero) {
				for (int j = 0; j < lookup.size(); j++) {
					newLookup[j]->push_back(lookup[j]->getAbundance(i), lookup[j]->getGroup());
				}
				//if there is a bin label use it otherwise make one
				string binLabel = "Otu";
				string sbinNumber = toString(i+1);
				if (sbinNumber.length() < snumBins.length()) { 
					int diff = snumBins.length() - sbinNumber.length();
					for (int h = 0; h < diff; h++) { binLabel += "0"; }
				}
				binLabel += sbinNumber; 
				if (i < m->currentSharedBinLabels.size()) {  binLabel = m->currentSharedBinLabels[i]; }
				
				newBinLabels.push_back(binLabel);
			}
		}
		
		for (int j = 0; j < lookup.size(); j++) {  delete lookup[j];  }
		
		lookup = newLookup;
		m->currentSharedBinLabels = newBinLabels;
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "SharedRAbundFloatVector", "eliminateZeroOTUS");
		exit(1);
	}
}
/***********************************************************************/

