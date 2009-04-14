/*
 *  list.cpp
 *  
 *
 *  Created by Pat Schloss on 8/8/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 */

using namespace std;


#include "sabundvector.hpp"
#include "rabundvector.hpp"
#include "ordervector.hpp"
#include "listvector.hpp"


/***********************************************************************/

ListVector::ListVector() : DataVector(), maxRank(0), numBins(0), numSeqs(0){};

/***********************************************************************/

ListVector::ListVector(int n):	DataVector(), data(n, "") , maxRank(0), numBins(0), numSeqs(0){};

/***********************************************************************/

ListVector::ListVector(string id, vector<string> lv) : DataVector(id), data(lv){
	try {
		for(int i=0;i<data.size();i++){
			if(data[i] != ""){
				int binSize = getNumNames(data[i]);
				numBins = i+1;
				if(binSize > maxRank)	{	maxRank = binSize;	}
				numSeqs += binSize;
			}
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ListVector class Function ListVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ListVector class function ListVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/**********************************************************************/

ListVector::ListVector(ifstream& f) : DataVector(), maxRank(0), numBins(0), numSeqs(0) {
	try {
		int hold;
		f >> label >> hold;
	
		data.assign(hold, "");
		string inputData = "";
	
		for(int i=0;i<hold;i++){
			f >> inputData;
			set(i, inputData);
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ListVector class Function ListVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ListVector class function ListVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/***********************************************************************/

void ListVector::set(int binNumber, string seqNames){
	try {
		int nNames_old = getNumNames(data[binNumber]);
		data[binNumber] = seqNames;
		int nNames_new = getNumNames(seqNames);
	
		if(nNames_old == 0)			{	numBins++;				}
		if(nNames_new == 0)			{	numBins--;				}
		if(nNames_new > maxRank)	{	maxRank = nNames_new;	}
	
		numSeqs += (nNames_new - nNames_old);
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ListVector class Function set. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ListVector class function set. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/***********************************************************************/

string ListVector::get(int index){
	return data[index];
}

/***********************************************************************/

void ListVector::push_back(string seqNames){
	try {
		data.push_back(seqNames);
		int nNames = getNumNames(seqNames);
	
		numBins++;
	
		if(nNames > maxRank)	{	maxRank = nNames;	}
	
		numSeqs += nNames;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ListVector class Function push_back. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ListVector class function push_back. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/***********************************************************************/

void ListVector::resize(int size){
	data.resize(size);		
}

/***********************************************************************/

int ListVector::size(){
	return data.size();
}
/***********************************************************************/

void ListVector::clear(){
	numBins = 0;
	maxRank = 0;
	numSeqs = 0;
	return data.clear();
	
}

/***********************************************************************/

void ListVector::print(ostream& output){
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
		cout << "Standard Error: " << e.what() << " has occurred in the ListVector class Function print. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ListVector class function print. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}


/***********************************************************************/

RAbundVector ListVector::getRAbundVector(){
	try {
		RAbundVector rav;
	
		for(int i=0;i<data.size();i++){
			int binSize = getNumNames(data[i]);
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
		cout << "Standard Error: " << e.what() << " has occurred in the ListVector class Function getRAbundVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ListVector class function getRAbundVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/***********************************************************************/

SAbundVector ListVector::getSAbundVector(){
	try {
		SAbundVector sav(maxRank+1);
	
		for(int i=0;i<data.size();i++){
			int binSize = getNumNames(data[i]);	
			sav.set(binSize, sav.get(binSize) + 1);	
		}
		sav.set(0, 0);
		sav.setLabel(label);
	
		return sav;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ListVector class Function getSAbundVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ListVector class function getSAbundVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/***********************************************************************/

OrderVector ListVector::getOrderVector(map<string,int>* orderMap = NULL){
	
	try {
		if(orderMap == NULL){
			OrderVector ov;
		
			for(int i=0;i<data.size();i++){
				int binSize = getNumNames(data[i]);		
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
				int length = listOTU.size();
				
				string seqName="";
			
				for(int j=0;j<length;j++){
				
					if(listOTU[j] != ','){
						seqName += listOTU[j];
					}
					else{
						if(orderMap->count(seqName) == 0){
							cerr << seqName << " not found, check *.names file\n";
							exit(1);
						}
					
						ov.set((*orderMap)[seqName], i);
						seqName = "";
					}						
				}
			
				if(orderMap->count(seqName) == 0){
					cerr << seqName << " not found, check *.names file\n";
					exit(1);
				}
				ov.set((*orderMap)[seqName], i);	
			}
		
			ov.setLabel(label);
			ov.getNumBins();
		
			return ov;		
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ListVector class Function getOrderVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ListVector class function getOrderVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/***********************************************************************/
