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
	m = MothurOut::getInstance();
	m->openInputFile(fName, fileHandle);
	filename = fName;
	m->saveNextLabel = "";
	
}
/***********************************************************************/

InputData::~InputData(){
	fileHandle.close();
	m->saveNextLabel = "";
}

/***********************************************************************/

InputData::InputData(string fName, string orderFileName, string f) : format(f){
	try {
		m = MothurOut::getInstance();
		ifstream ofHandle;
		m->openInputFile(orderFileName, ofHandle);
		string name;

		int count = 0;
	
		while(ofHandle){
			ofHandle >> name;
			orderMap[name] = count;
			count++;
			m->gobble(ofHandle);
		}
		ofHandle.close();
	
		m->openInputFile(fName, fileHandle);
	}
	catch(exception& e) {
		m->errorOut(e, "InputData", "InputData");
		exit(1);
	}
}
/***********************************************************************/

ListVector* InputData::getListVector(){
	try {
		if(!fileHandle.eof()){
			if(format == "list") {
				list = new ListVector(fileHandle);
			}else{ list = NULL;  }
					
			m->gobble(fileHandle);
			return list;
		}
		else{
			return NULL;
		}
	}
	catch(exception& e) {
		m->errorOut(e, "InputData", "getListVector");
		exit(1);
	}
}

/***********************************************************************/
ListVector* InputData::getListVector(string label){
	try {
		ifstream in;
		string  thisLabel;
		m->openInputFile(filename, in);
		
		if(in){

			if (format == "list") {
			
				while (in.eof() != true) {
					
					list = new ListVector(in);
					thisLabel = list->getLabel();
					
					//if you are at the last label
					if (thisLabel == label) {  break;  }
					//so you don't loose this memory
					else {	delete list;	}
					m->gobble(in);
				}
			}else{ list = NULL;  }
			
			in.close();
			return list;
		}
		else{
			return NULL;
		}
	}
	catch(exception& e) {
		m->errorOut(e, "InputData", "getListVector");
		exit(1);
	}
}
/***********************************************************************/
ListVector* InputData::getListVector(string label, bool resetFP){
	try {
		string  thisLabel;
		fileHandle.clear();
		fileHandle.seekg(0);
		
		if(fileHandle){

			if (format == "list") {
			
				while (fileHandle.eof() != true) {
					
					list = new ListVector(fileHandle); m->gobble(fileHandle);
					thisLabel = list->getLabel();
					
					//if you are at the last label
					if (thisLabel == label) {  break;  }
					//so you don't loose this memory
					else {	delete list;	}
				}
			}else{ list = NULL;  }
		
			return list;
		}
		else{
			return NULL;
		}
	}
	catch(exception& e) {
		m->errorOut(e, "InputData", "getListVector");
		exit(1);
	}
}

/***********************************************************************/

SharedListVector* InputData::getSharedListVector(){
	try {
		if(fileHandle){
			if (format == "shared")  {
				SharedList = new SharedListVector(fileHandle);
			}else{ SharedList = NULL;  }
					
			m->gobble(fileHandle);
			return SharedList;
		}
		else{
			return NULL;
		}
	}
	catch(exception& e) {
		m->errorOut(e, "InputData", "getSharedListVector");
		exit(1);
	}
}
/***********************************************************************/

SharedListVector* InputData::getSharedListVector(string label){
	try {
		ifstream in;
		string  thisLabel;
		m->openInputFile(filename, in);
		
		if(in){

			if (format == "shared")  {
			
				while (in.eof() != true) {
					
					SharedList = new SharedListVector(in);
					thisLabel = SharedList->getLabel();
					
					//if you are at the last label
					if (thisLabel == label) {  break;  }
					//so you don't loose this memory
					else {	delete SharedList;	}
					m->gobble(in);
				}

			}else{ SharedList = NULL;  }
				
			in.close();
			return SharedList;
			
		}else{
			return NULL;
		}
	}
	catch(exception& e) {
		m->errorOut(e, "InputData", "getSharedListVector");
		exit(1);
	}
}



/***********************************************************************/

SharedOrderVector* InputData::getSharedOrderVector(){
	try {
		if(fileHandle){
			if (format == "sharedfile")  {
				SharedOrder = new SharedOrderVector(fileHandle);
			}else{ SharedOrder = NULL;  }
				
			m->gobble(fileHandle);
			return SharedOrder;
			
		}else{
			return NULL;
		}
	}
	catch(exception& e) {
		m->errorOut(e, "InputData", "getSharedOrderVector");
		exit(1);
	}
}

/***********************************************************************/

SharedOrderVector* InputData::getSharedOrderVector(string label){
	try {
		ifstream in;
		string  thisLabel;
		m->openInputFile(filename, in);
		
		if(in){

			if (format == "sharedfile")  {
			
				while (in.eof() != true) {
					
					SharedOrder = new SharedOrderVector(in);
					thisLabel = SharedOrder->getLabel();
					
					//if you are at the last label
					if (thisLabel == label) {  break;  }
					//so you don't loose this memory
					else {	delete SharedOrder;	}
					m->gobble(in);
				}

			}else{ SharedOrder = NULL;  }
				
			in.close();
			return SharedOrder;
			
		}else{
			return NULL;
		}
	}
	catch(exception& e) {
		m->errorOut(e, "InputData", "getSharedOrderVector");
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
						
			m->gobble(fileHandle);
			
			output = new OrderVector();	
			*output = (input->getOrderVector());
		
			return output;
		}
		else{
			return NULL;
		}
	}
	catch(exception& e) {
		m->errorOut(e, "InputData", "getOrderVector");
		exit(1);
	}
}

/***********************************************************************/
OrderVector* InputData::getOrderVector(string label){
	try {
	
		ifstream in;
		string  thisLabel;
		m->openInputFile(filename, in);
		
		if(in){
			if((format == "list") || (format == "listorder")) {
			
				while (in.eof() != true) {
					
					input = new ListVector(in);
					thisLabel = input->getLabel();
					
					//if you are at the last label
					if (thisLabel == label) {  break;  }
					//so you don't loose this memory
					else {	delete input;	}
					m->gobble(in);
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
					m->gobble(in);
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
					m->gobble(in);
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
					m->gobble(in);
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
					m->gobble(in);
					
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
		m->errorOut(e, "InputData", "getOrderVector");
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
			m->gobble(fileHandle);
		}
				
		//this is created to signal to calling function that the input file is at eof
		vector<SharedRAbundVector*> null;  null.push_back(NULL);
		return null;
		
	}
	catch(exception& e) {
		m->errorOut(e, "InputData", "getSharedRAbundVectors");
		exit(1);
	}
}
/***********************************************************************/
vector<SharedRAbundVector*> InputData::getSharedRAbundVectors(string label){
	try {
		ifstream in;
		string  thisLabel;
		
		m->openInputFile(filename, in);
		m->saveNextLabel = "";
	
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
					m->gobble(in);
					
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
					m->gobble(in);
					
				}
			
			}
		}
				
		//this is created to signal to calling function that the input file is at eof
		vector<SharedRAbundVector*> null;  null.push_back(NULL);
		in.close();
		return null;
	
	}
	catch(exception& e) {
		m->errorOut(e, "InputData", "getSharedRAbundVectors");
		exit(1);
	}
}

/***********************************************************************/
//this is used when you don't need the order vector
vector<SharedRAbundFloatVector*> InputData::getSharedRAbundFloatVectors(){
	try {
		if(fileHandle){
			if (format == "relabund")  {
				SharedRAbundFloatVector* SharedRelAbund = new SharedRAbundFloatVector(fileHandle);
				if (SharedRelAbund != NULL) {
					return SharedRelAbund->getSharedRAbundFloatVectors();
				}
			}else if (format == "sharedfile")  {
				SharedRAbundVector* SharedRAbund = new SharedRAbundVector(fileHandle);
				if (SharedRAbund != NULL) {
					vector<SharedRAbundVector*> lookup = SharedRAbund->getSharedRAbundVectors(); 
					vector<SharedRAbundFloatVector*> lookupFloat = SharedRAbund->getSharedRAbundFloatVectors(lookup); 
					for (int i = 0; i < lookup.size(); i++) { delete lookup[i]; } lookup.clear();
					return lookupFloat;  
				}
						
			}
			m->gobble(fileHandle);
		}
				
		//this is created to signal to calling function that the input file is at eof
		vector<SharedRAbundFloatVector*> null;  null.push_back(NULL);
		return null;
		
	}
	catch(exception& e) {
		m->errorOut(e, "InputData", "getSharedRAbundFloatVectors");
		exit(1);
	}
}
/***********************************************************************/
vector<SharedRAbundFloatVector*> InputData::getSharedRAbundFloatVectors(string label){
	try {
		ifstream in;
		string  thisLabel;
		
		m->openInputFile(filename, in);
		m->saveNextLabel = "";
		
		if(in){
			if (format == "relabund")  {
				while (in.eof() != true) {
					
					SharedRAbundFloatVector* SharedRelAbund = new SharedRAbundFloatVector(in);
					if (SharedRelAbund != NULL) {
						thisLabel = SharedRelAbund->getLabel();
						//if you are at the last label
						if (thisLabel == label) {  in.close(); return SharedRelAbund->getSharedRAbundFloatVectors();  }
						else {
							//so you don't loose this memory
							vector<SharedRAbundFloatVector*> lookup = SharedRelAbund->getSharedRAbundFloatVectors(); 
							for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }
							delete SharedRelAbund;
						}
					}else{  break;  }
					m->gobble(in);
				}
			}else if (format == "sharedfile")  {
				while (in.eof() != true) {
					
					SharedRAbundVector* SharedRAbund = new SharedRAbundVector(in);
					if (SharedRAbund != NULL) {
						thisLabel = SharedRAbund->getLabel();
						
						//if you are at the last label
						if (thisLabel == label) {  
							in.close(); 
							vector<SharedRAbundVector*> lookup = SharedRAbund->getSharedRAbundVectors(); 
							vector<SharedRAbundFloatVector*> lookupFloat = SharedRAbund->getSharedRAbundFloatVectors(lookup); 
							for (int i = 0; i < lookup.size(); i++) { delete lookup[i]; } lookup.clear();
							return lookupFloat;  
						}else {
							//so you don't loose this memory
							vector<SharedRAbundVector*> lookup = SharedRAbund->getSharedRAbundVectors(); 
							for (int i = 0; i < lookup.size(); i++) { delete lookup[i]; } lookup.clear();
							delete SharedRAbund;
						}
					}else{  break;  }
					m->gobble(in);
				}
			}	
		}
		
				
		//this is created to signal to calling function that the input file is at eof
		vector<SharedRAbundFloatVector*> null;  null.push_back(NULL);
		in.close();
		return null;
	
	}
	catch(exception& e) {
		m->errorOut(e, "InputData", "getSharedRAbundFloatVectors");
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
			m->gobble(fileHandle);

			sabund = new SAbundVector();
			*sabund = (input->getSAbundVector());

			return sabund;
		}
		else{
			return NULL;
		}
	}
	catch(exception& e) {
		m->errorOut(e, "InputData", "getSAbundVector");
		exit(1);
	}
}
/***********************************************************************/
SAbundVector* InputData::getSAbundVector(string label){
	try {
	
		ifstream in;
		string  thisLabel;
		m->openInputFile(filename, in);
		
		if(in){
			if (format == "list") {
			
				while (in.eof() != true) {
					
					input = new ListVector(in);
					thisLabel = input->getLabel();
					
					//if you are at the last label
					if (thisLabel == label) {  break;  }
					//so you don't loose this memory
					else {	delete input;	}
					m->gobble(in);
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
					m->gobble(in);
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
					m->gobble(in);
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
					m->gobble(in);
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
					m->gobble(in);
					
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
		m->errorOut(e, "InputData", "getSAbundVector");
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
			
			m->gobble(fileHandle);

			rabund = new RAbundVector();
			*rabund = (input->getRAbundVector());

			return rabund;
		}
		else{
			return NULL;
		}
	}
	catch(exception& e) {
		m->errorOut(e, "InputData", "getRAbundVector");
		exit(1);
	}
}
/***********************************************************************/
RAbundVector* InputData::getRAbundVector(string label){
	try {
	
		ifstream in;
		string  thisLabel;
		m->openInputFile(filename, in);
		
		if(in){
			if (format == "list") {
			
				while (in.eof() != true) {
					
					input = new ListVector(in);
					thisLabel = input->getLabel();
					
					//if you are at the last label
					if (thisLabel == label) {  break;  }
					//so you don't loose this memory
					else {	delete input;	}
					m->gobble(in);
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
					m->gobble(in);
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
					m->gobble(in);
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
					m->gobble(in);
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
					m->gobble(in);
					
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
		m->errorOut(e, "InputData", "getRAbundVector");
		exit(1);
	}
}

/***********************************************************************/



