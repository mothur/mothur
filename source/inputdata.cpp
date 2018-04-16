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
			util.gobble(ofHandle);
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
                if (list != NULL) {
                    //pass labels to others distances in file
                    if (currentLabels.size() == 0) { currentLabels = list->getLabels(); }
                    else { list->setLabels(currentLabels);  }
                }
			}else{ list = NULL;  }
					
			util.gobble(fileHandle);
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
		util.openInputFile(filename, in);
        nextDistanceLabel = "";
		
		if(in){

			if (format == "list") {
			
				while (in.eof() != true) {
					
					list = new ListVector(in, nextDistanceLabel, otuTag);
					nextDistanceLabel = list->getLabel();
                    
                    if (list != NULL) {
                        //pass labels to others distances in file
                        if (currentLabels.size() == 0) { currentLabels = list->getLabels(); }
                        else { list->setLabels(currentLabels);  }
                    }
					
					//if you are at the last label
					if (nextDistanceLabel == label) {  break;  }
					//so you don't loose this memory
					else {	delete list;	}
					util.gobble(in);
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
		fileHandle.clear();
		fileHandle.seekg(0);
        nextDistanceLabel = "";
		
		if(fileHandle){

			if (format == "list") {
			
				while (fileHandle.eof() != true) {
					
					list = new ListVector(fileHandle, nextDistanceLabel, otuTag); util.gobble(fileHandle);
					nextDistanceLabel = list->getLabel();
                    
                     if (list != NULL) {
                         //pass labels to others distances in file
                         if (currentLabels.size() == 0) { currentLabels = list->getLabels(); }
                         else { list->setLabels(currentLabels);  }
                     }
					
					//if you are at the last label
					if (nextDistanceLabel == label) {  break;  }
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
				SharedList = new SharedListVector(fileHandle, groups, nextDistanceLabel, otuTag);
                if (SharedList != NULL) {
                    //pass labels to others distances in file
                    if (currentLabels.size() == 0) { currentLabels = SharedList->getLabels(); }
                    else { SharedList->setLabels(currentLabels);  }
                }

			}else{ SharedList = NULL;  }
					
			util.gobble(fileHandle);
			return SharedList;
		}
		else{ return NULL; }
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
		util.openInputFile(filename, in);
		
		if(in){

			if (format == "shared")  {
			
				while (!in.eof()) {
					
					SharedList = new SharedListVector(in, groups, nextDistanceLabel, otuTag);
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
					util.gobble(in);
				}

			}else{ SharedList = NULL;  }
				
			in.close();
			return SharedList;
			
		}else{ return NULL; }
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
                if (SharedOrder->getNumBins() == 0) { delete SharedOrder; SharedOrder =  NULL; } //no valid groups
			}else{ SharedOrder = NULL;  }
				
			util.gobble(fileHandle);
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
		util.openInputFile(filename, in);
        nextDistanceLabel = "";
		
		if(in){

			if (format == "sharedfile")  {
			
				while (in.eof() != true) {
					
					SharedOrder = new SharedOrderVector(in, groups, nextDistanceLabel);
					thisLabel = SharedOrder->getLabel();
                    if (SharedOrder->getNumBins() == 0) { delete SharedOrder; SharedOrder =  NULL; break; } //no valid groups
					
					//if you are at the last label
					if (thisLabel == label) {  break;  }
					//so you don't loose this memory
					else {	delete SharedOrder;	}
					util.gobble(in);
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
						
			util.gobble(fileHandle);
			
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
		util.openInputFile(filename, in);
		
		if(in){
			if((format == "list") || (format == "listorder")) {
                nextDistanceLabel = "";
				while (in.eof() != true) {
					
					input = new ListVector(in, nextDistanceLabel, otuTag);
					thisLabel = input->getLabel();
					
					//if you are at the last label
					if (thisLabel == label) {  break;  }
					//so you don't loose this memory
					else {	delete input;	}
					util.gobble(in);
				}
			}
			else if (format == "shared")  {
				nextDistanceLabel = "";
				while (in.eof() != true) {
					
					input = new SharedListVector(in, groups, nextDistanceLabel, otuTag);
					thisLabel = input->getLabel();
					
					//if you are at the last label
					if (thisLabel == label) {  break;  }
					//so you don't loose this memory
					else {	delete input;	}
					util.gobble(in);
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
					util.gobble(in);
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
					util.gobble(in);
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
					util.gobble(in);
					
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
                SharedRAbundVectors* shared = new SharedRAbundVectors(fileHandle, groups, nextDistanceLabel, otuTag);
                if (shared != NULL) {
                    //pass labels to others distances in file
                    if (currentLabels.size() == 0) { currentLabels = shared->getOTUNames(); }
                    else { shared->setOTUNames(currentLabels);  }
                    if (shared->getNumBins() == 0) { delete shared; shared = NULL; } //no valid groups
                }
                return shared;
            }else if (format == "shared") {
                SharedList = new SharedListVector(fileHandle, groups, nextDistanceLabel, otuTag);
                
                if (SharedList != NULL) {
                    //pass labels to others distances in file
                    if (currentLabels.size() == 0) { currentLabels = SharedList->getLabels(); }
                    else { SharedList->setLabels(currentLabels);  }

                    return SharedList->getSharedRAbundVector();
                }
            }
            util.gobble(fileHandle);
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
		
		util.openInputFile(filename, in);
		nextDistanceLabel = "";
	
		if(in){
			if (format == "sharedfile")  {
				while (in.eof() != true) {
					
					SharedRAbundVectors* SharedRAbund = new SharedRAbundVectors(in, groups, nextDistanceLabel, otuTag);
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
					util.gobble(in);
					
				}
			}else if (format == "shared") {
				while (in.eof() != true) {
					
					SharedList = new SharedListVector(in, groups, nextDistanceLabel, otuTag);
					
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
					util.gobble(in);
					
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
				SharedRAbundFloatVectors* SharedRelAbund = new SharedRAbundFloatVectors(fileHandle, groups, nextDistanceLabel, otuTag);
                if (SharedRelAbund->getNumBins() == 0) { delete SharedRelAbund; SharedRelAbund = NULL;  } //no valid groups
                else { //pass labels to others distances in file
                    if (currentLabels.size() == 0) { currentLabels = SharedRelAbund->getOTUNames(); }
                    else { SharedRelAbund->setOTUNames(currentLabels);  }
                }
                return SharedRelAbund;
			}else if (format == "sharedfile")  {
				SharedRAbundVectors* SharedRAbund = new SharedRAbundVectors(fileHandle, groups, nextDistanceLabel, otuTag);
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
			util.gobble(fileHandle);
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
		
		util.openInputFile(filename, in);
		nextDistanceLabel = "";
		
		if(in){
			if (format == "relabund")  {
				while (in.eof() != true) {
					
					SharedRAbundFloatVectors* SharedRelAbund = new SharedRAbundFloatVectors(in, groups, nextDistanceLabel, otuTag);
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
					util.gobble(in);
				}
			}else if (format == "sharedfile")  {
				while (in.eof() != true) {
					
					SharedRAbundVectors* SharedRAbund = new SharedRAbundVectors(in, groups, nextDistanceLabel, otuTag);
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
					util.gobble(in);
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
			util.gobble(fileHandle);

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
		util.openInputFile(filename, in);
		
		if(in){
			if (format == "list") {
                nextDistanceLabel = "";
				while (in.eof() != true) {
					
					input = new ListVector(in, nextDistanceLabel, otuTag);
					thisLabel = input->getLabel();
					
					//if you are at the last label
					if (thisLabel == label) {  break;  }
					//so you don't loose this memory
					else {	delete input;	}
					util.gobble(in);
				}
			}
			else if (format == "shared")  {
				nextDistanceLabel = "";
				while (in.eof() != true) {
					
					input = new SharedListVector(in, groups, nextDistanceLabel, otuTag);
					thisLabel = input->getLabel();
					
					//if you are at the last label
					if (thisLabel == label) {  break;  }
					//so you don't loose this memory
					else {	delete input;	}
					util.gobble(in);
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
					util.gobble(in);
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
					util.gobble(in);
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
					util.gobble(in);
					
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
			
			util.gobble(fileHandle);

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
		util.openInputFile(filename, in);
		
		if(in){
			if (format == "list") {
                nextDistanceLabel = "";
			
				while (in.eof() != true) {
					
					input = new ListVector(in, nextDistanceLabel, otuTag);
					thisLabel = input->getLabel();
					
					//if you are at the last label
					if (thisLabel == label) {  break;  }
					//so you don't loose this memory
					else {	delete input;	}
					util.gobble(in);
				}
			}
			else if (format == "shared")  {
                nextDistanceLabel = "";
				
				while (in.eof() != true) {
					
					input = new SharedListVector(in, groups, nextDistanceLabel, otuTag);
					thisLabel = input->getLabel();
					
					//if you are at the last label
					if (thisLabel == label) {  break;  }
					//so you don't loose this memory
					else {	delete input;	}
					util.gobble(in);
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
					util.gobble(in);
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
					util.gobble(in);
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
					util.gobble(in);
					
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



