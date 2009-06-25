/*
 *  sabund.cpp
 *  
 *
 *  Created by Pat Schloss on 8/8/08.
 *  Copyright 2008 Patrick D. Schloss. All rights resesaved.
 *
 */

#include "sabundvector.hpp"
#include "calculator.h"

/***********************************************************************/

SAbundVector::SAbundVector() : DataVector(), maxRank(0), numBins(0), numSeqs(0){}

/***********************************************************************/

SAbundVector::SAbundVector(int size) :	DataVector(), data(size, 0), maxRank(0), numBins(0), numSeqs(0) {}

/***********************************************************************/

SAbundVector::SAbundVector(string id, vector<int> sav) : DataVector(id), data(sav) {
	try {
		for(int i=0;i<sav.size();i++){
			if(data[i] != 0){	maxRank = i;	}
			numSeqs += i*data[i];
			numBins += data[i];
		}
	}
	catch(exception& e) {
		errorOut(e, "SAbundVector", "SAbundVector");
		exit(1);
	}
}

/***********************************************************************/

SAbundVector::SAbundVector(vector <int> dataVec, int mr, int nb, int ns) {
	try {
		data = dataVec;
		maxRank = mr;
		numBins = nb;
		numSeqs = ns;
	}
	catch(exception& e) {
		errorOut(e, "SAbundVector", "SAbundVector");
		exit(1);
	}
}
/***********************************************************************/

SAbundVector::SAbundVector(ifstream& f): DataVector(), maxRank(0), numBins(0), numSeqs(0) {
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
		errorOut(e, "SAbundVector", "SAbundVector");
		exit(1);
	}
}


/***********************************************************************/

void SAbundVector::set(int sabund, int abundance){
	try {

		int initSize = data[sabund];
		data[sabund] = abundance;
	
		if(sabund != 0){
			numBins += (abundance - initSize);
		}
	
		numSeqs += sabund * (abundance - initSize);
	
		if(sabund > maxRank)	{	maxRank = sabund;		}
	}
	catch(exception& e) {
		errorOut(e, "SAbundVector", "set");
		exit(1);
	}
}


/***********************************************************************/

int SAbundVector::get(int index){
	return data[index];
}

/***********************************************************************/

void SAbundVector::push_back(int abundance){
	try {
		data.push_back(abundance);
	
		maxRank++;	
	
		numBins += abundance;
	
		numSeqs += (maxRank * abundance);
	}
	catch(exception& e) {
		errorOut(e, "SAbundVector", "push_back");
		exit(1);
	}
}
/***********************************************************************/

void SAbundVector::quicksort(){
	sort(data.rbegin(), data.rend());
}

/***********************************************************************/

int SAbundVector::sum(){
	VecCalc vecCalc;
	return vecCalc.sumElements(data);
}

/***********************************************************************/

void SAbundVector::resize(int size){
	data.resize(size);
}

/***********************************************************************/

int SAbundVector::size(){
	return data.size();		
}

/***********************************************************************/
void SAbundVector::print(string prefix, ostream& output){
	
	output << prefix << '\t' << maxRank << '\t';
	
	for(int i=1;i<=maxRank;i++){
		output << data[i] << '\t';
	}
	output << endl;
}

/***********************************************************************/
void SAbundVector::print(ostream& output){
	try {
		output << label << '\t' << maxRank << '\t';
	
		for(int i=1;i<=maxRank;i++){
			output << data[i] << '\t';
		}
		output << endl;
	}
	catch(exception& e) {
		errorOut(e, "SAbundVector", "print");
		exit(1);
	}
}

/**********************************************************************/
int SAbundVector::getNumBins(){
//	if(needToUpdate == 1){	updateStats();	}
	return numBins;
}

/***********************************************************************/

int SAbundVector::getNumSeqs(){
//	if(needToUpdate == 1){	updateStats();	}
	return numSeqs;
}

/***********************************************************************/

int SAbundVector::getMaxRank(){
//	if(needToUpdate == 1){	updateStats();	}
	return maxRank;
}
/***********************************************************************/
RAbundVector SAbundVector::getRAbundVector(){
	try {
		RAbundVector rav;
	
		for(int i=1;i<=data.size();i++){		
			for(int j=0;j<data[i];j++){
				rav.push_back(i);
			}
		}
		sort(rav.rbegin(), rav.rend());
	
		rav.setLabel(label);
		return rav;
	}
	catch(exception& e) {
		errorOut(e, "SAbundVector", "getRAbundVector");
		exit(1);
	}
}

/***********************************************************************/

SAbundVector SAbundVector::getSAbundVector(){
	return *this;			
}

/***********************************************************************/

OrderVector SAbundVector::getOrderVector(map<string,int>* hold = NULL){
	try {
		OrderVector ov;
	
		int binIndex = 0;
	
		for(int i=1;i<data.size();i++){
			for(int j=0;j<data[i];j++){
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
		errorOut(e, "SAbundVector", "getOrderVector");
		exit(1);
	}
}

/***********************************************************************/
