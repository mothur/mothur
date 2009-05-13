/*
 *  inputdata.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 11/18/08.
 *  Copyright 2008 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "inputdata.h"
#include "ordervector.hpp"
#include "listvector.hpp"
#include "rabundvector.hpp"

/***********************************************************************/

InputData::InputData(string fName, string f) : format(f){
	
	openInputFile(fName, fileHandle);
	
}

/***********************************************************************/


InputData::~InputData(){
	
//	delete output;
	
}

/***********************************************************************/

InputData::InputData(string fName, string orderFileName, string f) : format(f){
	try {
		
		ifstream ofHandle;
		openInputFile(orderFileName, ofHandle);
		string name;

		int count = 0;
	
		while(ofHandle){
			ofHandle >> name;
			orderMap[name] = count;
			count++;
			gobble(ofHandle);
		}
		ofHandle.close();
	
		openInputFile(fName, fileHandle);
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the InputData class Function InputData. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the InputData class function InputData. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}
/***********************************************************************/

ListVector* InputData::getListVector(){
	try {
		if(fileHandle){
			if(format == "list") {
				list = new ListVector(fileHandle);
			}
					
			gobble(fileHandle);
			return list;
		}
		else{
			return 0;
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the InputData class Function getOrderVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the InputData class function getOrderVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

/***********************************************************************/

SharedListVector* InputData::getSharedListVector(){
	try {
		if(fileHandle){
			if (format == "shared")  {
				SharedList = new SharedListVector(fileHandle);
			}
					
			gobble(fileHandle);
			return SharedList;
		}
		else{
			return 0;
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the InputData class Function getOrderVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the InputData class function getOrderVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

/***********************************************************************/

SharedOrderVector* InputData::getSharedOrderVector(){
	try {
		if(fileHandle){
			if (format == "sharedfile")  {
				SharedOrder = new SharedOrderVector(fileHandle);
			}
				
			gobble(fileHandle);
			return SharedOrder;
			
		}else{
			return 0;
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the InputData class Function getSharedOrderVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the InputData class function getSharedOrderVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}



/***********************************************************************/

OrderVector* InputData::getOrderVector(){
	try {
		if(fileHandle){
			if(format == "list") {
				input = new ListVector(fileHandle);
			}
			else if (format == "shared")  {
				input = new SharedListVector(fileHandle);
			}
			else if(format == "rabund"){
				input = new RAbundVector(fileHandle);
			}
			else if(format == "order"){		
				input = new OrderVector(fileHandle);
			}
			else if(format == "sabund"){
				input = new SAbundVector(fileHandle);
			}
			else if(format == "listorder"){					
				input = new ListVector(fileHandle); 
			}
		
			gobble(fileHandle);
			output = new OrderVector();
			*output = (input->getOrderVector());
			
			return output;
		}
		else{
			return 0;
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the InputData class Function getOrderVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the InputData class function getOrderVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}
/***********************************************************************/
//this is used when you don't need the order vector
vector<SharedRAbundVector*> InputData::getSharedRAbundVectors(){
	try {
		if(fileHandle){
			if (format == "sharedfile")  {
				SharedRAbundVector* SharedRAbund = new SharedRAbundVector(fileHandle);
				if (SharedRAbund != NULL) {
					return SharedRAbund->getSharedRAbundVectors();
				}
			}else if (format == "shared") {
				SharedList = new SharedListVector(fileHandle);
				if (SharedList != NULL) {
					return SharedList->getSharedRAbundVector();
				}
			}
			gobble(fileHandle);
		}
				
		//this is created to signal to calling function that the input file is at eof
		vector<SharedRAbundVector*> null;  null.push_back(NULL);
		return null;
		
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the InputData class Function getSharedRAbundVectors. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the InputData class function getSharedRAbundVectors. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}


/***********************************************************************/

SAbundVector* InputData::getSAbundVector(){
	try {
		if(fileHandle){
			if (format == "list") {
				input = new ListVector(fileHandle);
			}
			else if (format == "shared")  {
				input = new SharedListVector(fileHandle);
			}
			else if(format == "rabund"){
				input = new RAbundVector(fileHandle);
			}
			else if(format == "order"){			
				input = new OrderVector(fileHandle);
			}
			else if(format == "sabund"){
				input = new SAbundVector(fileHandle);
			}
					
			gobble(fileHandle);

			sabund = new SAbundVector();
			*sabund = (input->getSAbundVector());

			return sabund;
		}
		else{
			return 0;
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the InputData class Function getSAbundVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the InputData class function getSAbundVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

/***********************************************************************/
RAbundVector* InputData::getRAbundVector(){
	try {
		if(fileHandle){
			if (format == "list") {
				input = new ListVector(fileHandle);
			}
			else if (format == "shared")  {
				input = new SharedListVector(fileHandle);
			}
			else if(format == "rabund"){
				input = new RAbundVector(fileHandle);
			}
			else if(format == "order"){			
				input = new OrderVector(fileHandle);
			}
			else if(format == "sabund"){
				input = new SAbundVector(fileHandle);
			}
					
			gobble(fileHandle);

			rabund = new RAbundVector();
			*rabund = (input->getRAbundVector());

			return rabund;
		}
		else{
			return NULL;
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the InputData class Function getRAbundVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the InputData class function getRAbundVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}
/***********************************************************************/



