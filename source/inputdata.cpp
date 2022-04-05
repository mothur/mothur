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
#include "sharedclrvectors.hpp"

/***********************************************************************/

InputData::InputData(string fName, string f, vector<string> userGroups) : format(f){
	m = MothurOut::getInstance();
	util.openInputFile(fName, fileHandle);
	filename = fName;
	nextDistanceLabel = "";
    groups = userGroups;
    otuTag = util.getTag(fName);
}
/***********************************************************************/

InputData::~InputData(){
	fileHandle.close();
	nextDistanceLabel = "";
}

/***********************************************************************/

InputData::InputData(string fName, string orderFileName, string f) : format(f){
	try {
		m = MothurOut::getInstance();
		ifstream ofHandle;
		util.openInputFile(orderFileName, ofHandle);
		string name;

		int count = 0;
	
		while(ofHandle){
			ofHandle >> name;
			orderMap[name] = count;
			count++;
			gobble(ofHandle);
		}
		ofHandle.close();
	
		util.openInputFile(fName, fileHandle);
		nextDistanceLabel = "";
        otuTag = util.getTag(fName);
		
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
				list = new ListVector(fileHandle, nextDistanceLabel, otuTag);
                if (list != nullptr) {
                    //pass labels to others distances in file
                    if (currentLabels.size() == 0) { currentLabels = list->getLabels(); }
                    else { list->setLabels(currentLabels);  }
                }
			}else{ list = nullptr;  }
					
			gobble(fileHandle);
			return list;
		}
		else{
			return nullptr;
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
		ifstream in; util.openInputFile(filename, in);
        nextDistanceLabel = "";
		
		if(in){

			if (format == "list") {
			
				while (!in.eof()) {
					
					list = new ListVector(in, nextDistanceLabel, otuTag);
					nextDistanceLabel = list->getLabel();
                    
                    if (list != nullptr) {
                        //pass labels to others distances in file
                        if (currentLabels.size() == 0) { currentLabels = list->getLabels(); }
                        else { list->setLabels(currentLabels);  }
                    }
					
					//if you are at the last label
					if (nextDistanceLabel == label) {  break;  }
					//so you don't loose this memory
					else {	delete list;	}
					gobble(in);
				}
			}else{ list = nullptr;  }
			
			in.close();
			return list;
		}
		else{
			return nullptr;
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
		fileHandle.clear();
		fileHandle.seekg(0);
        nextDistanceLabel = "";
		
		if(fileHandle){

			if (format == "list") {
			
				while (fileHandle.eof() != true) {
					
					list = new ListVector(fileHandle, nextDistanceLabel, otuTag); gobble(fileHandle);
					nextDistanceLabel = list->getLabel();
                    
                     if (list != nullptr) {
                         //pass labels to others distances in file
                         if (currentLabels.size() == 0) { currentLabels = list->getLabels(); }
                         else { list->setLabels(currentLabels);  }
                     }
					
					//if you are at the label you want
                    if (nextDistanceLabel == label) {  return list;  }
					else {	delete list;	} //so you don't loose this memory
				}
                
			}else{ return nullptr;   }
		}
        
        return nullptr;
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
				SharedList = new SharedListVector(fileHandle, groups, nextDistanceLabel, otuTag);
                if (SharedList != nullptr) {
                    //pass labels to others distances in file
                    if (currentLabels.size() == 0) { currentLabels = SharedList->getLabels(); }
                    else { SharedList->setLabels(currentLabels);  }
                }

			}else{ SharedList = nullptr;  }
					
			gobble(fileHandle);
			return SharedList;
		}
		else{ return nullptr; }
	}
	catch(exception& e) {
		m->errorOut(e, "InputData", "getSharedListVector");
		exit(1);
	}
}
/***********************************************************************/

SharedListVector* InputData::getSharedListVector(string label){
	try {
		
		string  thisLabel;
        ifstream in; util.openInputFile(filename, in);
        nextDistanceLabel = "";
		
		if(in){

			if (format == "shared")  {
			
				while (!in.eof()) {
					
					SharedList = new SharedListVector(in, groups, nextDistanceLabel, otuTag);
					thisLabel = SharedList->getLabel();
                    
                    if (SharedList != nullptr) {
                        //pass labels to others distances in file
                        if (currentLabels.size() == 0) { currentLabels = SharedList->getLabels(); }
                        else { SharedList->setLabels(currentLabels);  }
                    }
					//if you are at the last label
					if (thisLabel == label) {  break;  }
					//so you don't loose this memory
					else {	delete SharedList;	}
					gobble(in);
				}

			}else{ SharedList = nullptr;  }
				
			in.close();
			return SharedList;
			
		}else{ return nullptr; }
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
				SharedOrder = new SharedOrderVector(fileHandle, groups, nextDistanceLabel);
                if (SharedOrder->getNumBins() == 0) { delete SharedOrder; SharedOrder =  nullptr; } //no valid groups
			}else{ SharedOrder = nullptr;  }
				
			gobble(fileHandle);
			return SharedOrder;
			
		}else{
			return nullptr;
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
		string  thisLabel;
        ifstream in; util.openInputFile(filename, in);
        nextDistanceLabel = "";
		
		if(in){

			if (format == "sharedfile")  {
			
				while (!in.eof()) {
					
					SharedOrder = new SharedOrderVector(in, groups, nextDistanceLabel);
					thisLabel = SharedOrder->getLabel();
                    if (SharedOrder->getNumBins() == 0) { delete SharedOrder; SharedOrder =  nullptr; break; } //no valid groups
					
					//if you are at the last label
					if (thisLabel == label) {  break;  }
					//so you don't loose this memory
					else {	delete SharedOrder;	}
					gobble(in);
				}

			}else{ SharedOrder = nullptr;  }
				
			in.close();
			return SharedOrder;
			
		}else{
			return nullptr;
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
				input = new ListVector(fileHandle, nextDistanceLabel, otuTag);
			}
			else if (format == "shared")  {
				input = new SharedListVector(fileHandle, groups, nextDistanceLabel, otuTag);
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
			return nullptr;
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
		string  thisLabel;
        ifstream in; util.openInputFile(filename, in);
        nextDistanceLabel = "";
		
		if(in){
			if((format == "list") || (format == "listorder")) {
                nextDistanceLabel = "";
				while (!in.eof()) {
					
					input = new ListVector(in, nextDistanceLabel, otuTag);
					thisLabel = input->getLabel();
					
					//if you are at the last label
					if (thisLabel == label) {  break;  }
					//so you don't loose this memory
					else {	delete input;	}
					gobble(in);
				}
			}
			else if (format == "shared")  {
				nextDistanceLabel = "";
				while (!in.eof()) {
					
					input = new SharedListVector(in, groups, nextDistanceLabel, otuTag);
					thisLabel = input->getLabel();
					
					//if you are at the last label
					if (thisLabel == label) {  break;  }
					//so you don't loose this memory
					else {	delete input;	}
					gobble(in);
				}

			}
			else if(format == "rabund"){
				
				while (!in.eof()) {
					
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
				
				while (!in.eof()) {
					
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
				
				while (!in.eof()) {
					
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
			return nullptr;
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
                SharedRAbundVectors* shared = new SharedRAbundVectors(fileHandle, groups, nextDistanceLabel, otuTag);
                if (shared != nullptr) {
                    //pass labels to others distances in file
                    if (currentLabels.size() == 0) { currentLabels = shared->getOTUNames(); }
                    else { shared->setOTUNames(currentLabels);  }
                    if (shared->getNumBins() == 0) { delete shared; shared = nullptr; } //no valid groups
                }
                return shared;
            }else if (format == "shared") {
                SharedList = new SharedListVector(fileHandle, groups, nextDistanceLabel, otuTag);
                
                if (SharedList != nullptr) {
                    //pass labels to others distances in file
                    if (currentLabels.size() == 0) { currentLabels = SharedList->getLabels(); }
                    else { SharedList->setLabels(currentLabels);  }

                    return SharedList->getSharedRAbundVector();
                }
            }
            gobble(fileHandle);
        }
        
        //this is created to signal to calling function that the input file is at eof
        SharedRAbundVectors* null; null = nullptr;
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
		string  thisLabel;
		
        ifstream in; util.openInputFile(filename, in);
		nextDistanceLabel = "";
	
		if(in){
			if (format == "sharedfile")  {
				while (!in.eof()) {
					
					SharedRAbundVectors* SharedRAbund = new SharedRAbundVectors(in, groups, nextDistanceLabel, otuTag);
					if (SharedRAbund != nullptr) {
						thisLabel = SharedRAbund->getLabel();
                        
                        if (SharedRAbund->getNumBins() == 0) { delete SharedRAbund; SharedRAbund = nullptr; break; } //no valid groups
                        
                        //pass labels to others distances in file
                        if (currentLabels.size() == 0) { currentLabels = SharedRAbund->getOTUNames(); }
                        else { SharedRAbund->setOTUNames(currentLabels);  }
                        

						//if you are at the last label
						if (thisLabel == label) {  in.close(); return SharedRAbund;  }
						else {
							delete SharedRAbund;
						}
					}else{  break;  }
					gobble(in);
					
				}
			}else if (format == "shared") {
				while (!in.eof()) {
					
					SharedList = new SharedListVector(in, groups, nextDistanceLabel, otuTag);
					
					if (SharedList != nullptr) {
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
					gobble(in);
					
				}
			
			}
		}
				
		//this is created to signal to calling function that the input file is at eof
		SharedRAbundVectors* null;  null = (nullptr);
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
				SharedRAbundFloatVectors* SharedRelAbund = new SharedRAbundFloatVectors(fileHandle, groups, nextDistanceLabel, otuTag);
                if (SharedRelAbund->getNumBins() == 0) { delete SharedRelAbund; SharedRelAbund = nullptr;  } //no valid groups
                else { //pass labels to others distances in file
                    if (currentLabels.size() == 0) { currentLabels = SharedRelAbund->getOTUNames(); }
                    else { SharedRelAbund->setOTUNames(currentLabels);  }
                }
                return SharedRelAbund;
			}else if (format == "sharedfile")  {
				SharedRAbundVectors* SharedRAbund = new SharedRAbundVectors(fileHandle, groups, nextDistanceLabel, otuTag);
				if (SharedRAbund != nullptr) {
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
			gobble(fileHandle);
		}
				
		//this is created to signal to calling function that the input file is at eof
		SharedRAbundFloatVectors* null;  null = (nullptr);
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
		string  thisLabel;
		
        ifstream in; util.openInputFile(filename, in);
		nextDistanceLabel = "";
		
		if(in){
			if (format == "relabund")  {
				while (!in.eof()) {
					
					SharedRAbundFloatVectors* SharedRelAbund = new SharedRAbundFloatVectors(in, groups, nextDistanceLabel, otuTag);
					if (SharedRelAbund != nullptr) {
						thisLabel = SharedRelAbund->getLabel();
                        
                        if (SharedRelAbund->getNumBins() == 0) { delete SharedRelAbund; SharedRelAbund = nullptr; break;  } //no valid groups
                        
                        //pass labels to others distances in file
                        if (currentLabels.size() == 0) { currentLabels = SharedRelAbund->getOTUNames(); }
                        else { SharedRelAbund->setOTUNames(currentLabels);  }
						//if you are at the last label
						if (thisLabel == label) {  in.close(); return SharedRelAbund;  }
						else { delete SharedRelAbund; }
					}else{  break;  }
					gobble(in);
				}
			}else if (format == "sharedfile")  {
				while (!in.eof()) {
					
					SharedRAbundVectors* SharedRAbund = new SharedRAbundVectors(in, groups, nextDistanceLabel, otuTag);
					if (SharedRAbund != nullptr) {
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
					gobble(in);
				}
			}	
		}
		
				
		//this is created to signal to calling function that the input file is at eof
		SharedRAbundFloatVectors* null;  null = (nullptr);
		in.close();
		return null;
	
	}
	catch(exception& e) {
		m->errorOut(e, "InputData", "getSharedRAbundFloatVectors");
		exit(1);
	}
}
/***********************************************************************/
//this is used when you don't need the order vector
SharedCLRVectors* InputData::getSharedCLRVectors(){
    try {
        if(fileHandle){
            if (format == "clrfile")  {
                SharedCLRVectors* SharedCLR = new SharedCLRVectors(fileHandle, groups, nextDistanceLabel, otuTag);
                if (SharedCLR->getNumBins() == 0) { delete SharedCLR; SharedCLR = nullptr;  } //no valid groups
                else { //pass labels to others distances in file
                    if (currentLabels.size() == 0) { currentLabels = SharedCLR->getOTUNames(); }
                    else { SharedCLR->setOTUNames(currentLabels);  }
                }
                return SharedCLR;
            }
            gobble(fileHandle);
        }
                
        //this is created to signal to calling function that the input file is at eof
        SharedCLRVectors* null;  null = (nullptr);
        return null;
    }
    catch(exception& e) {
        m->errorOut(e, "InputData", "getSharedCLRVectors");
        exit(1);
    }
}
/***********************************************************************/
SharedCLRVectors* InputData::getSharedCLRVectors(string label){
    try {
        string  thisLabel;
        
        ifstream in; util.openInputFile(filename, in);
        nextDistanceLabel = "";
        
        if(in){
            if (format == "clrfile")  {
                while (!in.eof()) {
                    
                    SharedCLRVectors* SharedCLR = new SharedCLRVectors(in, groups, nextDistanceLabel, otuTag);
                    if (SharedCLR != nullptr) {
                        thisLabel = SharedCLR->getLabel();
                        
                        if (SharedCLR->getNumBins() == 0) { delete SharedCLR; SharedCLR = nullptr; break;  } //no valid groups
                        
                        //pass labels to others distances in file
                        if (currentLabels.size() == 0) { currentLabels = SharedCLR->getOTUNames(); }
                        else { SharedCLR->setOTUNames(currentLabels);  }
                        //if you are at the last label
                        if (thisLabel == label) {  in.close(); return SharedCLR;  }
                        else { delete SharedCLR; }
                    }else{  break;  }
                    gobble(in);
                }
            }
        }
        
        //this is created to signal to calling function that the input file is at eof
        SharedCLRVectors* null;  null = (nullptr);
        in.close();
        return null;
    
    }
    catch(exception& e) {
        m->errorOut(e, "InputData", "getSharedCLRVectors");
        exit(1);
    }
}
/***********************************************************************/

SAbundVector* InputData::getSAbundVector(){
	try {
		if(fileHandle){
			if (format == "list") {
				input = new ListVector(fileHandle, nextDistanceLabel, otuTag);
			}
			else if (format == "shared")  {
				input = new SharedListVector(fileHandle, groups, nextDistanceLabel, otuTag);
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
			return nullptr;
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
		string  thisLabel;
        ifstream in; util.openInputFile(filename, in);
        nextDistanceLabel = "";
		
		if(in){
			if (format == "list") {
                nextDistanceLabel = "";
				while (!in.eof()) {
					
					input = new ListVector(in, nextDistanceLabel, otuTag);
					thisLabel = input->getLabel();
					
					//if you are at the last label
					if (thisLabel == label) {  break;  }
					//so you don't loose this memory
					else {	delete input;	}
					gobble(in);
				}
			}
			else if (format == "shared")  {
				nextDistanceLabel = "";
				while (!in.eof()) {
					
					input = new SharedListVector(in, groups, nextDistanceLabel, otuTag);
					thisLabel = input->getLabel();
					
					//if you are at the last label
					if (thisLabel == label) {  break;  }
					//so you don't loose this memory
					else {	delete input;	}
					gobble(in);
				}

			}
			else if(format == "rabund"){
				
				while (!in.eof()) {
					
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
				
				while (!in.eof()) {
					
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
				
				while (!in.eof()) {
					
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
			return nullptr;
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
				input = new ListVector(fileHandle, nextDistanceLabel, otuTag);
			}
			else if (format == "shared")  {
				input = new SharedListVector(fileHandle, groups, nextDistanceLabel, otuTag);
			}
			else if(format == "rabund"){
				input = new RAbundVector(fileHandle);
			}
			else if(format == "order"){			
				input = new OrderVector(fileHandle);
			}
			else if(format == "sabund"){
				input = new SAbundVector(fileHandle);
            }else if (format == "sharedfile")  {
                SharedRAbundVectors* shared = new SharedRAbundVectors(fileHandle, groups, nextDistanceLabel, otuTag);
                if (shared != nullptr) {
                    //pass labels to others distances in file
                    if (currentLabels.size() == 0) { currentLabels = shared->getOTUNames(); }
                    else { shared->setOTUNames(currentLabels);  }
                    if (shared->getNumBins() == 0) { delete shared; shared = nullptr;  return nullptr; } //no valid groups
                }
                
                gobble(fileHandle);
                
                rabund = new RAbundVector();
                *rabund = (shared->getRAbundVector());
                
                delete shared;
                return rabund;
            }

			gobble(fileHandle);
            
			rabund = new RAbundVector();
			*rabund = (input->getRAbundVector());
            
            delete input;
			return rabund;
		}
		else{
			return nullptr;
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
		string  thisLabel;
        ifstream in; util.openInputFile(filename, in);
        nextDistanceLabel = "";
		
		if(in){
			if (format == "list") {
                nextDistanceLabel = "";
			
				while (!in.eof()) {
					
					input = new ListVector(in, nextDistanceLabel, otuTag);
					thisLabel = input->getLabel();
					
					//if you are at the last label
					if (thisLabel == label) {  break;  }
					//so you don't loose this memory
					else {	delete input;	}
					gobble(in);
				}
			}
			else if (format == "shared")  {
                nextDistanceLabel = "";
				
				while (!in.eof()) {
					
					input = new SharedListVector(in, groups, nextDistanceLabel, otuTag);
					thisLabel = input->getLabel();
					
					//if you are at the last label
					if (thisLabel == label) {  break;  }
					//so you don't loose this memory
					else {	delete input;	}
					gobble(in);
				}

			}
			else if(format == "rabund"){
				
				while (!in.eof()) {
					
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
				
				while (!in.eof()) {
					
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
				
				while (!in.eof()) {
					
					input = new SAbundVector(in);
					thisLabel = input->getLabel();
					
					//if you are at the last label
					if (thisLabel == label) {  break;  }
					//so you don't loose this memory
					else {	delete input;	}
					gobble(in);
					
				}

            }else if (format == "sharedfile")  {
                while (!in.eof()) {
                    
                    SharedRAbundVectors* shared = new SharedRAbundVectors(in, groups, nextDistanceLabel, otuTag);
                    if (shared != nullptr) {
                        thisLabel = shared->getLabel();
                        
                        if (shared->getNumBins() == 0) { delete shared; shared = nullptr; in.close(); return nullptr; } //no valid groups
                        
                        //pass labels to others distances in file
                        if (currentLabels.size() == 0) { currentLabels = shared->getOTUNames(); }
                        else { shared->setOTUNames(currentLabels);  }
                        
                        
                        //if you are at the last label
                        if (thisLabel == label) {
                            in.close();
                            
                            rabund = new RAbundVector();
                            *rabund = (shared->getRAbundVector());
                            
                            delete shared;
                            return rabund;
                        }
                        else { delete shared;  }
                    }else{ in.close(); return nullptr;  }
                    gobble(in);
                }
                
            }
			
			
			in.close();		

			rabund = new RAbundVector();
			*rabund = (input->getRAbundVector());

			return rabund;
		}
		else{
			return nullptr;
		}
	}
	catch(exception& e) {
		m->errorOut(e, "InputData", "getRAbundVector");
		exit(1);
	}
}

/***********************************************************************/



