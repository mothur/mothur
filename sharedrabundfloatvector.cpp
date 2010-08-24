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

SharedRAbundFloatVector::SharedRAbundFloatVector() : DataVector(), maxRank(0.0), numBins(0), numSeqs(0.0) {globaldata = GlobalData::getInstance();}
/***********************************************************************/

SharedRAbundFloatVector::~SharedRAbundFloatVector() {}

/***********************************************************************/
SharedRAbundFloatVector::SharedRAbundFloatVector(int n) : DataVector(), maxRank(0.0), numBins(n), numSeqs(0.0) {
		globaldata = GlobalData::getInstance();
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
		globaldata = GlobalData::getInstance();
		
		if (globaldata->gGroupmap == NULL) {  groupmap = new GroupMap(); }
		
		int num, count;
		float inputData;
		count = 0;  
		string holdLabel, nextLabel, groupN;
		individualFloat newguy;
		
		for (int i = 0; i < lookup.size(); i++) {  delete lookup[i]; lookup[i] = NULL; }
		lookup.clear();
		
		//read in first row since you know there is at least 1 group.
		f >> label >> groupN >> num;
		holdLabel = label;
		
		//add new vector to lookup
		SharedRAbundFloatVector* temp = new SharedRAbundFloatVector();
		lookup.push_back(temp);
		lookup[0]->setLabel(label);
		lookup[0]->setGroup(groupN);
		
		if (globaldata->gGroupmap == NULL) { 
			//save group in groupmap
			groupmap->namesOfGroups.push_back(groupN);
			groupmap->groupIndex[groupN] = 0;
		}
		
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
			
			if (globaldata->gGroupmap == NULL) { 
				//save group in groupmap
				groupmap->namesOfGroups.push_back(groupN);
				groupmap->groupIndex[groupN] = count;
			}
			
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
		
		//put file pointer back since you are now at a new distance label
		for (int i = 0; i < nextLabel.length(); i++) { f.unget();  }
	
		if (globaldata->gGroupmap == NULL) { globaldata->gGroupmap = groupmap;  }
		
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
		m->errorOut(e, "SharedRAbundVector", "set");
		exit(1);
	}
}

/***********************************************************************/
float SharedRAbundFloatVector::getAbundance(int index){
	return data[index].abundance;	
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
void SharedRAbundFloatVector::print(ostream& output){
	try {
		output << numBins << '\t';
	
		for(int i=0;i<data.size();i++){		output << data[i].abundance << '\t';		}
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
/***********************************************************************/
vector<SharedRAbundFloatVector*> SharedRAbundFloatVector::getSharedRAbundFloatVectors(){
	try {
		SharedUtil* util;
		util = new SharedUtil();
		
		util->setGroups(globaldata->Groups, globaldata->gGroupmap->namesOfGroups);

		for (int i = 0; i < lookup.size(); i++) {
			//if this sharedrabund is not from a group the user wants then delete it.
			if (util->isValidGroup(lookup[i]->getGroup(), globaldata->Groups) == false) { 
				delete lookup[i]; lookup[i] = NULL;
				lookup.erase(lookup.begin()+i); 
				i--; 
			}
		}
		
		delete util;
	
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
/***********************************************************************/

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
/***********************************************************************/
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

/***********************************************************************/

