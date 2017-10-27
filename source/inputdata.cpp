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
#include "sharedrabundvectors.hpp"

/***********************************************************************/

InputData::InputData(string fName, string f, vector<string> userGroups) : format(f){
	m = MothurOut::getInstance();
	m->openInputFile(fName, fileHandle);
	filename = fName;
	m->setSaveNextLabel("");
    groups = userGroups;
}
/***********************************************************************/

InputData::~InputData(){
	fileHandle.close();
	m->setSaveNextLabel("");
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
		m->setSaveNextLabel("");
		
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
                if (list != NULL) {
                    //pass labels to others distances in file
                    if (currentLabels.size() == 0) { currentLabels = list->getLabels(); }
                    else { list->setLabels(currentLabels);  }
                }
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
        m->setSaveNextLabel("");
		
		if(in){

			if (format == "list") {
			
				while (in.eof() != true) {
					
					list = new ListVector(in);
					thisLabel = list->getLabel();
                    
                    if (list != NULL) {
                        //pass labels to others distances in file
                        if (currentLabels.size() == 0) { currentLabels = list->getLabels(); }
                        else { list->setLabels(currentLabels);  }
                    }
					
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
        m->setSaveNextLabel("");
		
		if(fileHandle){

			if (format == "list") {
			
				while (fileHandle.eof() != true) {
					
					list = new ListVector(fileHandle); m->gobble(fileHandle);
					thisLabel = list->getLabel();
                    
                     if (list != NULL) {
                         //pass labels to others distances in file
                         if (currentLabels.size() == 0) { currentLabels = list->getLabels(); }
                         else { list->setLabels(currentLabels);  }
                     }
					
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
				SharedList = new SharedListVector(fileHandle, groups);
                if (SharedList != NULL) {
                    //pass labels to others distances in file
                    if (currentLabels.size() == 0) { currentLabels = SharedList->getLabels(); }
                    else { SharedList->setLabels(currentLabels);  }
                }

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
					
					SharedList = new SharedListVector(in, groups);
					thisLabel = SharedList->getLabel();
                    
                    if (SharedList != NULL) {
                        //pass labels to others distances in file
                        if (currentLabels.size() == 0) { currentLabels = SharedList->getLabels(); }
                        else { SharedList->setLabels(currentLabels);  }
                    }
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
				SharedOrder = new SharedOrderVector(fileHandle, groups);
                if (SharedOrder->getNumBins() == 0) { delete SharedOrder; SharedOrder =  NULL; } //no valid groups
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
        m->setSaveNextLabel("");
		
		if(in){

			if (format == "sharedfile")  {
			
				while (in.eof() != true) {
					
					SharedOrder = new SharedOrderVector(in, groups);
					thisLabel = SharedOrder->getLabel();
                    if (SharedOrder->getNumBins() == 0) { delete SharedOrder; SharedOrder =  NULL; break; } //no valid groups
					
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
				input = new SharedListVector(fileHandle, groups);
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
                m->setSaveNextLabel("");
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
				m->setSaveNextLabel("");
				while (in.eof() != true) {
					
					input = new SharedListVector(in, groups);
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
SharedRAbundVectors* InputData::getSharedRAbundVectors(){
    try {
        if(fileHandle){
            if (format == "sharedfile")  {
                SharedRAbundVectors* shared = new SharedRAbundVectors(fileHandle, groups);
                if (shared != NULL) {
                    //pass labels to others distances in file
                    if (currentLabels.size() == 0) { currentLabels = shared->getOTUNames(); }
                    else { shared->setOTUNames(currentLabels);  }
                    if (shared->getNumBins() == 0) { delete shared; shared = NULL; } //no valid groups
                }
                return shared;
            }else if (format == "shared") {
                SharedList = new SharedListVector(fileHandle, groups);
                
                if (SharedList != NULL) {
                    //pass labels to others distances in file
                    if (currentLabels.size() == 0) { currentLabels = SharedList->getLabels(); }
                    else { SharedList->setLabels(currentLabels);  }

                    return SharedList->getSharedRAbundVector();
                }
            }
            m->gobble(fileHandle);
        }
        
        //this is created to signal to calling function that the input file is at eof
        SharedRAbundVectors* null; null = NULL;
        return null;
        
    }
    catch(exception& e) {
        m->errorOut(e, "InputData", "getSharedRAbundVectors");
        exit(1);
    }
}

/***********************************************************************/
SharedRAbundVectors* InputData::getSharedRAbundVectors(string label){
	try {
		ifstream in;
		string  thisLabel;
		
		m->openInputFile(filename, in);
		m->setSaveNextLabel("");
	
		if(in){
			if (format == "sharedfile")  {
				while (in.eof() != true) {
					
					SharedRAbundVectors* SharedRAbund = new SharedRAbundVectors(in, groups);
					if (SharedRAbund != NULL) {
						thisLabel = SharedRAbund->getLabel();
                        
                        if (SharedRAbund->getNumBins() == 0) { delete SharedRAbund; SharedRAbund = NULL; break; } //no valid groups
                        
                        //pass labels to others distances in file
                        if (currentLabels.size() == 0) { currentLabels = SharedRAbund->getOTUNames(); }
                        else { SharedRAbund->setOTUNames(currentLabels);  }
                        

						//if you are at the last label
						if (thisLabel == label) {  in.close(); return SharedRAbund;  }
						else {
							delete SharedRAbund;
						}
					}else{  break;  }
					m->gobble(in);
					
				}
			}else if (format == "shared") {
				while (in.eof() != true) {
					
					SharedList = new SharedListVector(in, groups);
					
					if (SharedList != NULL) {
						thisLabel = SharedList->getLabel();
                        //pass labels to others distances in file
                        if (currentLabels.size() == 0) { currentLabels = SharedList->getLabels(); }
                        else { SharedList->setLabels(currentLabels);  }

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
		SharedRAbundVectors* null;  null = (NULL);
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
SharedRAbundFloatVectors* InputData::getSharedRAbundFloatVectors(){
	try {
		if(fileHandle){
			if (format == "relabund")  {
				SharedRAbundFloatVectors* SharedRelAbund = new SharedRAbundFloatVectors(fileHandle, groups);
                if (SharedRelAbund->getNumBins() == 0) { delete SharedRelAbund; SharedRelAbund = NULL;  } //no valid groups
                else { //pass labels to others distances in file
                    if (currentLabels.size() == 0) { currentLabels = SharedRelAbund->getOTUNames(); }
                    else { SharedRelAbund->setOTUNames(currentLabels);  }
                }
                return SharedRelAbund;
			}else if (format == "sharedfile")  {
				SharedRAbundVectors* SharedRAbund = new SharedRAbundVectors(fileHandle, groups);
				if (SharedRAbund != NULL) {
                    //pass labels to others distances in file
                    if (currentLabels.size() == 0) { currentLabels = SharedRAbund->getOTUNames(); }
                    else { SharedRAbund->setOTUNames(currentLabels);  }

					vector<SharedRAbundFloatVector*> lookup = SharedRAbund->getSharedRAbundFloatVectors();
                    SharedRAbundFloatVectors* SharedRelAbund = new SharedRAbundFloatVectors();
                    SharedRelAbund->setOTUNames(currentLabels);
                    for (int i = 0; i < lookup.size(); i++) { SharedRelAbund->push_back(lookup[i]); }
					delete SharedRAbund;
                    return SharedRelAbund;
				}
						
			}
			m->gobble(fileHandle);
		}
				
		//this is created to signal to calling function that the input file is at eof
		SharedRAbundFloatVectors* null;  null = (NULL);
		return null;
	}
	catch(exception& e) {
		m->errorOut(e, "InputData", "getSharedRAbundFloatVectors");
		exit(1);
	}
}
/***********************************************************************/
SharedRAbundFloatVectors* InputData::getSharedRAbundFloatVectors(string label){
	try {
		ifstream in;
		string  thisLabel;
		
		m->openInputFile(filename, in);
		m->setSaveNextLabel("");
		
		if(in){
			if (format == "relabund")  {
				while (in.eof() != true) {
					
					SharedRAbundFloatVectors* SharedRelAbund = new SharedRAbundFloatVectors(in, groups);
					if (SharedRelAbund != NULL) {
						thisLabel = SharedRelAbund->getLabel();
                        
                        if (SharedRelAbund->getNumBins() == 0) { delete SharedRelAbund; SharedRelAbund = NULL; break;  } //no valid groups
                        
                        //pass labels to others distances in file
                        if (currentLabels.size() == 0) { currentLabels = SharedRelAbund->getOTUNames(); }
                        else { SharedRelAbund->setOTUNames(currentLabels);  }
						//if you are at the last label
						if (thisLabel == label) {  in.close(); return SharedRelAbund;  }
						else { delete SharedRelAbund; }
					}else{  break;  }
					m->gobble(in);
				}
			}else if (format == "sharedfile")  {
				while (in.eof() != true) {
					
					SharedRAbundVectors* SharedRAbund = new SharedRAbundVectors(in, groups);
					if (SharedRAbund != NULL) {
						thisLabel = SharedRAbund->getLabel();
                        
                        //pass labels to others distances in file
                        if (currentLabels.size() == 0) { currentLabels = SharedRAbund->getOTUNames(); }
                        else { SharedRAbund->setOTUNames(currentLabels);  }
						
						//if you are at the last label
						if (thisLabel == label) {  
							in.close(); 
                            vector<SharedRAbundFloatVector*> lookup = SharedRAbund->getSharedRAbundFloatVectors();
                            SharedRAbundFloatVectors* SharedRelAbund = new SharedRAbundFloatVectors();
                            SharedRelAbund->setOTUNames(currentLabels);
                            for (int i = 0; i < lookup.size(); i++) { SharedRelAbund->push_back(lookup[i]); }
                            delete SharedRAbund;
                            return SharedRelAbund;
						}else { delete SharedRAbund; }
					}else{  break;  }
					m->gobble(in);
				}
			}	
		}
		
				
		//this is created to signal to calling function that the input file is at eof
		SharedRAbundFloatVectors* null;  null = (NULL);
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
				input = new SharedListVector(fileHandle, groups);
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
                m->setSaveNextLabel("");
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
				m->setSaveNextLabel("");
				while (in.eof() != true) {
					
					input = new SharedListVector(in, groups);
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
				input = new SharedListVector(fileHandle, groups);
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
                m->setSaveNextLabel("");
			
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
                m->setSaveNextLabel("");
				
				while (in.eof() != true) {
					
					input = new SharedListVector(in, groups);
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



