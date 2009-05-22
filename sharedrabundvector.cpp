/*
 *  sharedvector.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 12/5/08.
 *  Copyright 2008 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */


using namespace std;

#include "sharedrabundvector.h" 
#include "sabundvector.hpp"
#include "ordervector.hpp"
#include "sharedutilities.h"


/***********************************************************************/

SharedRAbundVector::SharedRAbundVector() : DataVector(), maxRank(0), numBins(0), numSeqs(0) {};

/***********************************************************************/

SharedRAbundVector::SharedRAbundVector(int n) : DataVector(), maxRank(0), numBins(n), numSeqs(0) {
		individual newGuy;
		//initialize data
		for (int i=0; i< n; i++) {
			newGuy.bin = i;
			newGuy.abundance = 0;
			data.push_back(newGuy);
		}
};

/***********************************************************************

SharedRAbundVector::SharedRAbundVector(string id, vector<individual> rav) : DataVector(id), data(rav) {
	try {
		numBins = 0;
		maxRank = 0;
		numSeqs = 0;
		
		for(int i=0;i<data.size();i++){
			if(data[i].abundance != 0)		{	numBins = i+1;		}
			if(data[i].abundance > maxRank)	{	maxRank = data[i].abundance;	}
			numSeqs += data[i].abundance;
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SharedRAbundVector class Function SharedRAbundVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SharedRAbundVector class function SharedRAbundVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}


/***********************************************************************/
//reads a shared file
SharedRAbundVector::SharedRAbundVector(ifstream& f) : DataVector(), maxRank(0), numBins(0), numSeqs(0) {
	try {
		globaldata = GlobalData::getInstance();
		
		if (globaldata->gGroupmap == NULL) {  groupmap = new GroupMap(); }
		
		int num, inputData, pos, count;
		count = 0;  
		string holdLabel, nextLabel, groupN;
		individual newguy;
		
		for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }
		lookup.clear();
		
		//read in first row since you know there is at least 1 group.
		f >> label >> groupN >> num;
		holdLabel = label;
		
		//add new vector to lookup
		SharedRAbundVector* temp = new SharedRAbundVector();
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
			
			lookup[0]->push_back(inputData, i, groupN); //abundance, bin, group
			push_back(inputData, i, groupN);
			numSeqs += inputData;
			numBins++;
			if (inputData > maxRank) { maxRank = inputData; }
			
		}
		
		//save position in file in case next line is a new label.
		pos = f.tellg();
		
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
			temp = new SharedRAbundVector();
			lookup.push_back(temp);
			lookup[count]->setLabel(label);
			lookup[count]->setGroup(groupN);

			//fill vector.  
			for(int i=0;i<num;i++){
				f >> inputData;
				lookup[count]->push_back(inputData, i, groupN); //abundance, bin, group
			}
			
			//save position in file in case next line is a new label.
			pos = f.tellg();
	
			if (f.eof() != true) { f >> nextLabel; }
		}
		
		//put file pointer back since you are now at a new distance label
		f.seekg(pos, ios::beg);
	
		if (globaldata->gGroupmap == NULL) { globaldata->gGroupmap = groupmap; }
		
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SharedRAbundVector class Function SharedRAbundVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SharedRAbundVector class function SharedRAbundVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/***********************************************************************/

SharedRAbundVector::~SharedRAbundVector() {

}

/***********************************************************************/

void SharedRAbundVector::set(int binNumber, int newBinSize, string groupname){
	try {
		int oldBinSize = data[binNumber].abundance;
		data[binNumber].abundance = newBinSize;
		data[binNumber].group = groupname;
	
		if(newBinSize > maxRank)	{	maxRank = newBinSize;	}
	
		numSeqs += (newBinSize - oldBinSize);
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SharedRAbundVector class Function set. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SharedRAbundVector class function set. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
/***********************************************************************/

void SharedRAbundVector::setData(vector <individual> newData){
	data = newData;
}

/***********************************************************************/

int SharedRAbundVector::getAbundance(int index){
	return data[index].abundance;
	
}
/***********************************************************************/

int SharedRAbundVector::numNZ(){
	int sum = 0;
	for(int i = 1; i < numBins; i++)
		if(data[i].abundance > 0)
			sum++;
	return sum;
}
/***********************************************************************/

void SharedRAbundVector::sortD(){
	struct individual indObj;
	sort(data.begin()+1, data.end(), indObj);
}
/***********************************************************************/

individual SharedRAbundVector::get(int index){
	return data[index];
	
}
/***********************************************************************/

vector <individual> SharedRAbundVector::getData(){
	return data;
}
/***********************************************************************/

void SharedRAbundVector::push_back(int binSize, int otu, string groupName){
	try {
		individual newGuy;
		newGuy.abundance = binSize;
		newGuy.group = groupName;
		newGuy.bin = otu;
		
		data.push_back(newGuy);
		numBins++;
	
		if(binSize > maxRank){
			maxRank = binSize;
		}
	
		numSeqs += binSize;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SharedRAbundVector class Function push_back. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SharedRAbundVector class function push_back. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/***********************************************************************/

void SharedRAbundVector::insert(int binSize, int otu, string groupName){
	try {
		individual newGuy;
		newGuy.abundance = binSize;
		newGuy.group = groupName;
		newGuy.bin = otu;
		
		data.insert(data.begin()+otu, newGuy);
		numBins++;
	
		if(binSize > maxRank){
			maxRank = binSize;
		}
	
		numSeqs += binSize;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SharedRAbundVector class Function insert. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SharedRAbundVector class function insert. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/***********************************************************************/

void SharedRAbundVector::push_front(int binSize, int otu, string groupName){
	try {
		individual newGuy;
		newGuy.abundance = binSize;
		newGuy.group = groupName;
		newGuy.bin = otu;
		
		data.insert(data.begin(), newGuy);
		numBins++;
	
		if(binSize > maxRank){
			maxRank = binSize;
		}
	
		numSeqs += binSize;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SharedRAbundVector class Function push_front. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SharedRAbundVector class function push_front. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/***********************************************************************/
void SharedRAbundVector::pop_back(){
	numSeqs -= data[data.size()-1].abundance;
	numBins--;
	return data.pop_back();
}

/***********************************************************************/


vector<individual>::reverse_iterator SharedRAbundVector::rbegin(){
	return data.rbegin();				
}

/***********************************************************************/

vector<individual>::reverse_iterator SharedRAbundVector::rend(){
	return data.rend();					
}

/***********************************************************************/
void SharedRAbundVector::resize(int size){
	
	data.resize(size);
}

/***********************************************************************/

int SharedRAbundVector::size(){
	return data.size();
}

/***********************************************************************/
void SharedRAbundVector::print(ostream& output){
	try {
		output << numBins << '\t';
	
		for(int i=0;i<data.size();i++){		output << data[i].abundance << '\t';		}
		output << endl;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SharedRAbundVector class Function print. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SharedRAbundVector class function print. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
/***********************************************************************/
string SharedRAbundVector::getGroup(){
	return group;
}

/***********************************************************************/

void SharedRAbundVector::setGroup(string groupName){
	group = groupName;
}
/***********************************************************************/
int SharedRAbundVector::getGroupIndex()  { return index; }
/***********************************************************************/
void SharedRAbundVector::setGroupIndex(int vIndex)	{ index = vIndex; }
/***********************************************************************/
int SharedRAbundVector::getNumBins(){
	return numBins;
}

/***********************************************************************/

int SharedRAbundVector::getNumSeqs(){
	return numSeqs;
}

/***********************************************************************/

int SharedRAbundVector::getMaxRank(){
	return maxRank;
}
/***********************************************************************/

SharedRAbundVector SharedRAbundVector::getSharedRAbundVector(){
	return *this;			
}
/***********************************************************************/
vector<SharedRAbundVector*> SharedRAbundVector::getSharedRAbundVectors(){
	try {
		SharedUtil* util;
		util = new SharedUtil();
		
		util->setGroups(globaldata->Groups, globaldata->gGroupmap->namesOfGroups);

		for (int i = 0; i < lookup.size(); i++) {
			//if this sharedrabund is not from a group the user wants then delete it.
			if (util->isValidGroup(lookup[i]->getGroup(), globaldata->Groups) == false) { 
				delete lookup[i]; 
				lookup.erase(lookup.begin()+i); 
				i--; 
			}
		}

		return lookup;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SharedRAbundVector class Function getSharedRAbundVectors. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SharedRAbundVector class function getSharedRAbundVectors. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
/***********************************************************************/

RAbundVector SharedRAbundVector::getRAbundVector() {
	try {
		RAbundVector rav(data.size());
		
		for (int i = 0; i < data.size(); i++) {
			rav.set(i, data[i].abundance);
		}
	
		return rav;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SharedRAbundVector class Function getRAbundVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SharedRAbundVector class function getRAbundVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
/***********************************************************************/

RAbundVector SharedRAbundVector::getRAbundVector2() {
	try {
		RAbundVector rav;
		for(int i = 0; i < numBins; i++)
			if(data[i].abundance != 0)
				rav.push_back(data[i].abundance-1);
		return rav;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SharedRAbundVector class Function getRAbundVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SharedRAbundVector class function getRAbundVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
/***********************************************************************/

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
		cout << "Standard Error: " << e.what() << " has occurred in the SharedRAbundVector class Function getSharedSAbundVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SharedRAbundVector class function getSharedSAbundVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
/***********************************************************************/

SAbundVector SharedRAbundVector::getSAbundVector() {
	try {
		SAbundVector sav(maxRank+1);
		
		for(int i=0;i<data.size();i++){
			int abund = data[i].abundance;
			sav.set(abund, sav.get(abund) + 1);
		}
		sav.set(0, 0);
		sav.setLabel(label);
		return sav;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SharedRAbundVector class Function getSAbundVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SharedRAbundVector class function getSAbundVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/***********************************************************************/

SharedOrderVector SharedRAbundVector::getSharedOrderVector() {
	try {
		SharedOrderVector ov;
	
		for(int i=0;i<data.size();i++){
			for(int j=0;j<data[i].abundance;j++){
				ov.push_back(data[i].bin, data[i].abundance, data[i].group);
			}
		}
		random_shuffle(ov.begin(), ov.end());

		ov.setLabel(label);	
		ov.updateStats();
		
		return ov;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SharedRAbundVector class Function getOrderVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SharedRAbundVector class function getOrderVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
/***********************************************************************/

OrderVector SharedRAbundVector::getOrderVector(map<string,int>* nameMap = NULL) {
	try {
		OrderVector ov;
	
		for(int i=0;i<data.size();i++){
			for(int j=0;j<data[i].abundance;j++){
				ov.push_back(i);
			}
		}
		random_shuffle(ov.begin(), ov.end());

		ov.setLabel(label);	
		return ov;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SharedRAbundVector class Function getOrderVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SharedRAbundVector class function getOrderVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/***********************************************************************/

