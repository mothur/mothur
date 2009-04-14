/*
 *  order.cpp
 *  
 *
 *  Created by Pat Schloss on 8/8/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 */

using namespace std;

#include "ordervector.hpp"


/***********************************************************************/

OrderVector::OrderVector() : DataVector() {}

/***********************************************************************/

//OrderVector::OrderVector(int ns) : DataVector(), data(ns, -1) {};

/***********************************************************************/

OrderVector::OrderVector(string id, vector<int> ov) : 
											DataVector(id), data(ov)
{
	updateStats();	
}

/***********************************************************************/

OrderVector::OrderVector(ifstream& f) : DataVector() {
	try {
		int hold;
	
		f >> label;
		f >> hold;
	
		data.assign(hold, -1);
	
		int inputData;
	
		for(int i=0;i<hold;i++){
			f >> inputData;
			set(i, inputData);
		}
	
		updateStats();
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the OrderVector class Function OrderVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the OrderVector class function OrderVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/***********************************************************************/


int OrderVector::getNumBins(){
	if(needToUpdate == 1){	updateStats();	}
	return numBins;
}

/***********************************************************************/

int OrderVector::getNumSeqs(){
	if(needToUpdate == 1){	updateStats();	}
	return numSeqs;
}

/***********************************************************************/

int OrderVector::getMaxRank(){
	if(needToUpdate == 1){	updateStats();	}
	return maxRank;
}

/***********************************************************************/



void OrderVector::set(int index, int binNumber){
	
	data[index] = binNumber;
	needToUpdate = 1;
	
}

/***********************************************************************/

int OrderVector::get(int index){
	return data[index];			
}

/***********************************************************************/

void OrderVector::push_back(int index){

	data.push_back(index);
	needToUpdate = 1;
	
}

/***********************************************************************/

void OrderVector::print(ostream& output){
	try {
		output << label << '\t' << numSeqs << '\t';
	
		for(int i=0;i<data.size();i++){
			output << data[i] << '\t';
		}
		output << endl;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the OrderVector class Function print. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the OrderVector class function print. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/***********************************************************************/

void OrderVector::print(string prefix, ostream& output){
	try {
		output << prefix << '\t' << numSeqs << '\t';
	
		for(int i=0;i<numSeqs;i++){
			output << data[i] << '\t';
		}
		output << endl;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the OrderVector class Function print. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the OrderVector class function print. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

/***********************************************************************/

void OrderVector::resize(int){
	cout << "resize() did nothing in class OrderVector";
}

/***********************************************************************/

int OrderVector::size(){
	return data.size();					
}

/***********************************************************************/

vector<int>::iterator OrderVector::begin(){
	return data.begin();	
}

/***********************************************************************/

vector<int>::iterator OrderVector::end(){
	return data.end();		
}

/***********************************************************************/

RAbundVector OrderVector::getRAbundVector(){
	try {
		RAbundVector rav(data.size());
	
		for(int i=0;i<numSeqs;i++){
			rav.set(data[i], rav.get(data[i]) + 1);
		}	
		sort(rav.rbegin(), rav.rend());
		for(int i=numSeqs-1;i>=0;i--){
			if(rav.get(i) == 0){	rav.pop_back();	}
			else{
				break;
			}
		}
		rav.setLabel(label);

		return rav;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the OrderVector class Function getRAbundVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the OrderVector class function getRAbundVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

/***********************************************************************/

SAbundVector OrderVector::getSAbundVector(){
	
	RAbundVector rav(this->getRAbundVector());
	return rav.getSAbundVector();

}

/***********************************************************************/

OrderVector OrderVector::getOrderVector(map<string,int>* hold = 0){
	return *this;			
}

/***********************************************************************/

void OrderVector::updateStats(){
	try {
		needToUpdate = 0;
	//	int maxBinVectorLength = 0;
		numSeqs = 0;
		numBins = 0;
		maxRank = 0;
	
		for(int i=0;i<data.size();i++){
			if(data[i] != -1){
				numSeqs++;
			}
		}
	
		vector<int> hold(numSeqs);
	
		for(int i=0;i<numSeqs;i++){
			if(data[i] != -1){
				hold[data[i]] = hold[data[i]]+1;
			}
		}	
		for(int i=0;i<numSeqs;i++){
			if(hold[i] > 0)			{	numBins++;			}
			if(hold[i] > maxRank)	{	maxRank = hold[i];	}
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the OrderVector class Function updateStats. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the OrderVector class function updateStats. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

/***********************************************************************/

