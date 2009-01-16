/*
 *  sharedSharedSAbundVector.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 12/10/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "sharedsabundvector.h"
#include "sabundvector.hpp"
#include "datavector.hpp"
#include "utilities.hpp"
#include <exception>

using namespace std;


/***********************************************************************/

SharedSAbundVector::SharedSAbundVector() : DataVector(), maxRank(0), numBins(0), numSeqs(0){};

/***********************************************************************/

SharedSAbundVector::SharedSAbundVector(int size) :	DataVector(), maxRank(0), numBins(0), numSeqs(0) {
		individual newGuy;
		//initialize data
		for (int i=0; i< size; i++) {
			newGuy.bin = i;
			newGuy.abundance = 0;
			data.push_back(newGuy);
		}
	//	for(int i=0;i<data.size();i++){
	//		if(data[i].abundance != 0){	maxRank = i;	}
	//		numSeqs += i*data[i].abundance;
	//		numBins += data[i].abundance;
	//	}

};

/***********************************************************************

SharedSAbundVector::SharedSAbundVector(string id, vector<int> sav) : DataVector(id), data(sav) {
	try {
		
		for(int i=0;i<sav.size();i++){
			if(data[i] != 0){	maxRank = i;	}
			numSeqs += i*data[i];
			numBins += data[i];
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SharedSAbundVector class Function SharedSAbundVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SharedSAbundVector class function SharedSAbundVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/***********************************************************************

SharedSAbundVector::SharedSAbundVector(ifstream& f): DataVector(), maxRank(0), numBins(0), numSeqs(0) {
	try {
		int hold;
		f >> label >> hold;
	
		data.assign(hold+1, 0);
		int inputData;
	
		for(int i=1;i<=hold;i++){
			f >> inputData;
			set(i, inputData);
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SharedSAbundVector class Function SharedSAbundVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SharedSAbundVector class function SharedSAbundVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/***********************************************************************/

void SharedSAbundVector::set(int bin, int abundance, string groupName){
	try {

		int initSize = data[bin].abundance;
		data[bin].abundance = abundance;
		data[bin].group = groupName;
	
		if(bin != 0){
			numBins += (abundance - initSize);
		}
	
		numSeqs += bin * (abundance - initSize);
	
		if(bin > maxRank)	{	maxRank = bin;		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SharedSAbundVector class Function set. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SharedSAbundVector class function set. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

/***********************************************************************/

individual SharedSAbundVector::get(int index){
	return data[index];
}
/***********************************************************************/

int SharedSAbundVector::getAbundance(int index){
	return data[index].abundance;
}

/***********************************************************************/

void SharedSAbundVector::push_back(int abundance, int bin, string groupName){
	try {
		individual newGuy;
		newGuy.abundance = abundance;
		newGuy.bin = bin;
		newGuy.group = groupName;
		
		data.push_back(newGuy);
	
		maxRank++;	
	
		numBins += abundance;
	
		numSeqs += (maxRank * abundance);
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SharedSAbundVector class Function push_back. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SharedSAbundVector class function push_back. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		
}

/***********************************************************************/

void SharedSAbundVector::resize(int size){
	data.resize(size);
}

/***********************************************************************/

int SharedSAbundVector::size(){
	return data.size();		
}

/***********************************************************************/
void SharedSAbundVector::print(ostream& output){
	try {
		output << label << '\t' << maxRank << '\t';
	
		for(int i=1;i<=maxRank;i++){
			output << data[i].abundance << '\t';
		}
		output << endl;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SharedSAbundVector class Function print. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SharedSAbundVector class function print. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		
}
/***********************************************************************/
string SharedSAbundVector::getGroup(){
	return group;
}

/***********************************************************************/

void SharedSAbundVector::setGroup(string groupName){
	group = groupName;
}

/**********************************************************************/
int SharedSAbundVector::getNumBins(){
	return numBins;
}

/***********************************************************************/

int SharedSAbundVector::getNumSeqs(){
	return numSeqs;
}

/***********************************************************************/

int SharedSAbundVector::getMaxRank(){
	return maxRank;
}
/***********************************************************************/
RAbundVector SharedSAbundVector::getRAbundVector(){
	try {
		RAbundVector rav;
	
		for(int i=1;i<=data.size();i++){		
			for(int j=0;j<data[i].abundance;j++){
				rav.push_back(i);
			}
		}
		sort(rav.rbegin(), rav.rend());
	
		rav.setLabel(label);
		return rav;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SharedSAbundVector class Function getRAbundVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SharedSAbundVector class function getRAbundVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		
}
/***********************************************************************/
SAbundVector SharedSAbundVector::getSAbundVector(){
	try {
		RAbundVector rav;
		SAbundVector sav;
		
		rav = getRAbundVector();
		sav = rav.getSAbundVector();
		return sav;
	
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SharedSAbundVector class Function getSAbundVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SharedSAbundVector class function getSAbundVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		
}

/***********************************************************************/

bool compareMembers (individual member, individual member2){

  if(member.abundance < member2.abundance){
    return true;   }   
  else{
	return false; 
  }
}

/***********************************************************************/
SharedRAbundVector SharedSAbundVector::getSharedVector(){
	try {
		SharedRAbundVector rav;
		
		int binNumber = 0;
		for(int i=1;i<=data.size();i++){		
			for(int j=0;j<data[i].abundance;j++){
				rav.push_back(i, binNumber, data[i].group);
				binNumber++;
			}
		}
		sort(rav.rbegin(), rav.rend(), compareMembers);
	
		rav.setLabel(label);
		rav.setGroup(group);
		
		return rav;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SharedSAbundVector class Function getSharedVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SharedSAbundVector class function getSharedVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		
}


/***********************************************************************/

SharedSAbundVector SharedSAbundVector::getSharedSAbundVector(){
	return *this;			
}

/***********************************************************************/

OrderVector SharedSAbundVector::getOrderVector(map<string,int>* hold = NULL){
	try {
		OrderVector ov;
	
		int binIndex = 0;
	
		for(int i=1;i<data.size();i++){
			for(int j=0;j<data[i].abundance;j++){
				for(int k=0;k<i;k++){
					ov.push_back(binIndex);
				}
				binIndex++;
			}
		}
	
		random_shuffle(ov.begin(), ov.end());

		ov.setLabel(label);
		ov.getNumBins();
		return ov;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SharedSAbundVector class Function getOrderVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SharedSAbundVector class function getOrderVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		
}

/***********************************************************************/

