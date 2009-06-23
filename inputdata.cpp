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
	filename = fName;
	
}

/***********************************************************************/


InputData::~InputData(){
	fileHandle.close();
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
			return NULL;
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the InputData class Function getListVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the InputData class function getListVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

/***********************************************************************/

ListVector* InputData::getListVector(string label){
	try {
		ifstream in;
		string  thisLabel;
		openInputFile(filename, in);
		
		if(in){

			if (format == "list") {
			
				while (in.eof() != true) {
					
					list = new ListVector(in);
					thisLabel = list->getLabel();
					
					//if you are at the last label
					if (thisLabel == label) {  break;  }
					//so you don't loose this memory
					else {	delete list;	}
					gobble(in);
				}
			}
			
			in.close();
			return list;
		}
		else{
			return NULL;
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the InputData class Function getListVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the InputData class function getListVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
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
			return NULL;
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

SharedListVector* InputData::getSharedListVector(string label){
	try {
		ifstream in;
		string  thisLabel;
		openInputFile(filename, in);
		
		if(in){

			if (format == "shared")  {
			
				while (in.eof() != true) {
					
					SharedList = new SharedListVector(in);
					thisLabel = SharedList->getLabel();
					
					//if you are at the last label
					if (thisLabel == label) {  break;  }
					//so you don't loose this memory
					else {	delete SharedList;	}
					gobble(in);
				}

			}
				
			in.close();
			return SharedList;
			
		}else{
			return NULL;
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the InputData class Function getSharedListVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the InputData class function getSharedListVector. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
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
			return NULL;
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

SharedOrderVector* InputData::getSharedOrderVector(string label){
	try {
		ifstream in;
		string  thisLabel;
		openInputFile(filename, in);
		
		if(in){

			if (format == "sharedfile")  {
			
				while (in.eof() != true) {
					
					SharedOrder = new SharedOrderVector(in);
					thisLabel = SharedOrder->getLabel();
					
					//if you are at the last label
					if (thisLabel == label) {  break;  }
					//so you don't loose this memory
					else {	delete SharedOrder;	}
					gobble(in);
				}

			}
				
			in.close();
			return SharedOrder;
			
		}else{
			return NULL;
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
			if((format == "list") || (format == "listorder")) {
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
			output = new OrderVector();
			*output = (input->getOrderVector());
			
			return output;
		}
		else{
			return NULL;
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
OrderVector* InputData::getOrderVector(string label){
	try {
	
		ifstream in;
		string  thisLabel;
		openInputFile(filename, in);
		
		if(in){
			if((format == "list") || (format == "listorder")) {
			
				while (in.eof() != true) {
					
					input = new ListVector(in);
					thisLabel = input->getLabel();
					
					//if you are at the last label
					if (thisLabel == label) {  break;  }
					//so you don't loose this memory
					else {	delete input;	}
					gobble(in);
				}
			}
			else if (format == "shared")  {
				
				while (in.eof() != true) {
					
					input = new SharedListVector(in);
					thisLabel = input->getLabel();
					
					//if you are at the last label
					if (thisLabel == label) {  break;  }
					//so you don't loose this memory
					else {	delete input;	}
					gobble(in);
				}

			}
			else if(format == "rabund"){
				
				while (in.eof() != true) {
					
					input = new RAbundVector(in);
					thisLabel = input->getLabel();
					
					//if you are at the last label
					if (thisLabel == label) {  break;  }
					//so you don't loose this memory
					else {	delete input;	}
					gobble(in);
				}

			}
			else if(format == "order"){			
				
				while (in.eof() != true) {
					
					input = new OrderVector(in);
					thisLabel = input->getLabel();
					
					//if you are at the last label
					if (thisLabel == label) {  break;  }
					//so you don't loose this memory
					else {	delete input;	}
					gobble(in);
				}

			}
			else if(format == "sabund"){
				
				while (in.eof() != true) {
					
					input = new SAbundVector(in);
					thisLabel = input->getLabel();
					
					//if you are at the last label
					if (thisLabel == label) {  break;  }
					//so you don't loose this memory
					else {	delete input;	}
					gobble(in);
					
				}

			}
			
			in.close();		

			output = new OrderVector();
			*output = (input->getOrderVector());
			
			return output;

		}
		else{
			return NULL;
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
vector<SharedRAbundVector*> InputData::getSharedRAbundVectors(string label){
	try {
		ifstream in;
		string  thisLabel;
		
		openInputFile(filename, in);
		
		if(in){
			if (format == "sharedfile")  {
				while (in.eof() != true) {
					
					SharedRAbundVector* SharedRAbund = new SharedRAbundVector(in);
					if (SharedRAbund != NULL) {
						thisLabel = SharedRAbund->getLabel();
						//if you are at the last label
						if (thisLabel == label) {  in.close(); return SharedRAbund->getSharedRAbundVectors();  }
						else {
							//so you don't loose this memory
							vector<SharedRAbundVector*> lookup = SharedRAbund->getSharedRAbundVectors(); 
							for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }
							delete SharedRAbund;
						}
					}else{  break;  }
					gobble(in);
					
				}
			}else if (format == "shared") {
				while (in.eof() != true) {
					
					SharedList = new SharedListVector(in);
					if (SharedList != NULL) {
						thisLabel = SharedList->getLabel();
						//if you are at the last label
						if (thisLabel == label) {  in.close(); return SharedList->getSharedRAbundVector();  }
						else {
							//so you don't loose this memory
							delete SharedList;
						}
					}else{  break;  }
					gobble(in);
					
				}
			
			}
		}
				
		//this is created to signal to calling function that the input file is at eof
		vector<SharedRAbundVector*> null;  null.push_back(NULL);
		in.close();
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
			return NULL;
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
SAbundVector* InputData::getSAbundVector(string label){
	try {
	
		ifstream in;
		string  thisLabel;
		openInputFile(filename, in);
		
		if(in){
			if (format == "list") {
			
				while (in.eof() != true) {
					
					input = new ListVector(in);
					thisLabel = input->getLabel();
					
					//if you are at the last label
					if (thisLabel == label) {  break;  }
					//so you don't loose this memory
					else {	delete input;	}
					gobble(in);
				}
			}
			else if (format == "shared")  {
				
				while (in.eof() != true) {
					
					input = new SharedListVector(in);
					thisLabel = input->getLabel();
					
					//if you are at the last label
					if (thisLabel == label) {  break;  }
					//so you don't loose this memory
					else {	delete input;	}
					gobble(in);
				}

			}
			else if(format == "rabund"){
				
				while (in.eof() != true) {
					
					input = new RAbundVector(in);
					thisLabel = input->getLabel();
					
					//if you are at the last label
					if (thisLabel == label) {  break;  }
					//so you don't loose this memory
					else {	delete input;	}
					gobble(in);
				}

			}
			else if(format == "order"){			
				
				while (in.eof() != true) {
					
					input = new OrderVector(in);
					thisLabel = input->getLabel();
					
					//if you are at the last label
					if (thisLabel == label) {  break;  }
					//so you don't loose this memory
					else {	delete input;	}
					gobble(in);
				}

			}
			else if(format == "sabund"){
				
				while (in.eof() != true) {
					
					input = new SAbundVector(in);
					thisLabel = input->getLabel();
					
					//if you are at the last label
					if (thisLabel == label) {  break;  }
					//so you don't loose this memory
					else {	delete input;	}
					gobble(in);
					
				}

			}
			
			in.close();		

			sabund = new SAbundVector();
			*sabund = (input->getSAbundVector());

			return sabund;

		}
		else{
			return NULL;
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
RAbundVector* InputData::getRAbundVector(string label){
	try {
	
		ifstream in;
		string  thisLabel;
		openInputFile(filename, in);
		
		if(in){
			if (format == "list") {
			
				while (in.eof() != true) {
					
					input = new ListVector(in);
					thisLabel = input->getLabel();
					
					//if you are at the last label
					if (thisLabel == label) {  break;  }
					//so you don't loose this memory
					else {	delete input;	}
					gobble(in);
				}
			}
			else if (format == "shared")  {
				
				while (in.eof() != true) {
					
					input = new SharedListVector(in);
					thisLabel = input->getLabel();
					
					//if you are at the last label
					if (thisLabel == label) {  break;  }
					//so you don't loose this memory
					else {	delete input;	}
					gobble(in);
				}

			}
			else if(format == "rabund"){
				
				while (in.eof() != true) {
					
					input = new RAbundVector(in);
					thisLabel = input->getLabel();
					
					//if you are at the last label
					if (thisLabel == label) {  break;  }
					//so you don't loose this memory
					else {	delete input;	}
					gobble(in);
				}

			}
			else if(format == "order"){			
				
				while (in.eof() != true) {
					
					input = new OrderVector(in);
					thisLabel = input->getLabel();
					
					//if you are at the last label
					if (thisLabel == label) {  break;  }
					//so you don't loose this memory
					else {	delete input;	}
					gobble(in);
				}

			}
			else if(format == "sabund"){
				
				while (in.eof() != true) {
					
					input = new SAbundVector(in);
					thisLabel = input->getLabel();
					
					//if you are at the last label
					if (thisLabel == label) {  break;  }
					//so you don't loose this memory
					else {	delete input;	}
					gobble(in);
					
				}

			}
			
			in.close();		

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



