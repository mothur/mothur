//
//  oligos.cpp
//  Mothur
//
//  Created by Sarah Westcott on 4/4/14.
//  Copyright (c) 2014 Schloss Lab. All rights reserved.
//

#include "oligos.h"
#include "utils.hpp"

/**************************************************************************************************/

Oligos::Oligos(string o){
	try {
		m = MothurOut::getInstance();
        hasPPrimers = false; hasPBarcodes = false; pairedOligos = false; reversePairs = true;
        indexBarcode = 0; indexPairedBarcode = 0; indexPrimer = 0; indexPairedPrimer = 0;
		oligosfile = o;
        reversePairs = true;
        readOligos();
		if (pairedOligos) {
            numBarcodes = pairedBarcodes.size();
            numFPrimers = pairedPrimers.size();
        }else {
            numBarcodes = barcodes.size();
            numFPrimers = primers.size();
        }
	}
	catch(exception& e) {
		m->errorOut(e, "Oligos", "Oligos");
		exit(1);
	}
}
/**************************************************************************************************/

Oligos::Oligos(){
	try {
		m = MothurOut::getInstance();
        hasPPrimers = false; hasPBarcodes = false; pairedOligos = false; reversePairs = true;
        indexBarcode = 0; indexPairedBarcode = 0; indexPrimer = 0; indexPairedPrimer = 0;
        numFPrimers = 0; numBarcodes = 0;
	}
	catch(exception& e) {
		m->errorOut(e, "Oligos", "Oligos");
		exit(1);
	}
}
/**************************************************************************************************/
int Oligos::read(string o){
	try {
		oligosfile = o;
        readOligos();
        if (pairedOligos) {
            numBarcodes = pairedBarcodes.size();
            numFPrimers = pairedPrimers.size();
        }else {
            numBarcodes = barcodes.size();
            numFPrimers = primers.size();
        }
        return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "Oligos", "read");
		exit(1);
	}
}
/**************************************************************************************************/
int Oligos::read(string o, bool reverse){
	try {
		oligosfile = o;
        reversePairs = reverse;
        readOligos();
        if (pairedOligos) {
            numBarcodes = pairedBarcodes.size();
            numFPrimers = pairedPrimers.size();
        }else {
            numBarcodes = barcodes.size();
            numFPrimers = primers.size();
        }
        return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "Oligos", "read");
		exit(1);
	}
}
/**************************************************************************************************/
vector<string> Oligos::getSRAGroupNames(){
    try {
        vector<string> sraGroupNames;
        set<string> uniqueNames;
        
        if (pairedOligos) {
            for(map<int, oligosPair>::iterator itBar = pairedBarcodes.begin();itBar != pairedBarcodes.end();itBar++){
                for(map<int, oligosPair>::iterator itPrimer = pairedPrimers.begin();itPrimer != pairedPrimers.end(); itPrimer++){
                    
                    if (m->getControl_pressed()) { return sraGroupNames; }
                    
                    string primerName = getPrimerName(itPrimer->first);
                    string barcodeName = getBarcodeName(itBar->first);
                    
                    if ((primerName == "ignore") || (barcodeName == "ignore")) { } //do nothing
                    else if ((primerName == "") && (barcodeName == "")) { } //do nothing
                    else {
                        string comboGroupName = "";
                        string comboName = "";
                        
                        if(primerName == ""){
                            comboGroupName = barcodeName;
                        }else{
                            if(barcodeName == ""){
                                comboGroupName = primerName;
                            }
                            else{
                                comboGroupName = barcodeName + "." + primerName;
                            }
                        }
                        
                        if(((itPrimer->second).forward+(itPrimer->second).reverse) == ""){
                            if ((itBar->second).forward != "NONE") { comboName += (itBar->second).forward; }
                            if ((itBar->second).reverse != "NONE") {
                                if (comboName == "") {  comboName += (itBar->second).reverse; }
                                else {  comboName += ("."+(itBar->second).reverse);  }
                            }
                        }else{
                            if(((itBar->second).forward+(itBar->second).reverse) == ""){
                                if ((itPrimer->second).forward != "NONE") { comboName += (itPrimer->second).forward; }
                                if ((itPrimer->second).reverse != "NONE") {
                                    if (comboName == "") {  comboName += (itPrimer->second).reverse; }
                                    else {  comboName += ("."+(itPrimer->second).reverse);  }
                                }
                            }
                            else{
                                if ((itBar->second).forward != "NONE") { comboName += (itBar->second).forward; }
                                if ((itBar->second).reverse != "NONE") {
                                    if (comboName == "") {  comboName += (itBar->second).reverse; }
                                    else {  comboName += ("."+(itBar->second).reverse);  }
                                }
                                if ((itPrimer->second).forward != "NONE") {
                                    if (comboName == "") {  comboName += (itPrimer->second).forward; }
                                    else {  comboName += ("."+(itPrimer->second).forward);  }
                                }
                                if ((itPrimer->second).reverse != "NONE") {
                                    if (comboName == "") {  comboName += (itPrimer->second).reverse; }
                                    else {  comboName += ("."+(itPrimer->second).reverse);  }
                                }
                            }
                        }
                        
                        if (comboName != "") {  comboGroupName +=  "_" + comboName;  }
                        uniqueNames.insert(comboGroupName);
                    }
                }
            }
        }else {
            
            for(map<string, int>::iterator itBar = barcodes.begin();itBar != barcodes.end();itBar++){
                for(map<string, int>::iterator itPrimer = primers.begin();itPrimer != primers.end(); itPrimer++){
                    
                    string primerName = getPrimerName(itPrimer->second);
                    string barcodeName = getBarcodeName(itBar->second);
                    
                    if ((primerName == "ignore") || (barcodeName == "ignore")) { } //do nothing
                    else if ((primerName == "") && (barcodeName == "")) { } //do nothing
                    else {
                        string comboGroupName = "";
                        string comboName = "";
                        
                        if(primerName == ""){
                            comboGroupName = barcodeName;
                        }else{
                            if(barcodeName == ""){
                                comboGroupName = primerName;
                            }
                            else{
                                comboGroupName = barcodeName + "." + primerName;
                            }
                        }
                        
                        if(itPrimer->first == ""){
                            comboName = itBar->first;
                        }else{
                            if(itBar->first == ""){
                                comboName = itPrimer->first;
                            }
                            else{
                                comboName = itBar->first + "." + itPrimer->first;
                            }
                        }
                        
                        if (comboName != "") {  comboGroupName +=  "_" + comboName;  }
                        uniqueNames.insert(comboGroupName);
                    }
                }
            }
        }
        
        if (uniqueNames.size() == 0) {
            m->mothurOut("[ERROR]: your oligos file does not contain any group names.\n");  m->setControl_pressed(true);
        }else {
            if (m->getDebug()) { int count = 0; for (set<string>::iterator it = uniqueNames.begin(); it != uniqueNames.end(); it++) { m->mothurOut("[DEBUG]: " + toString(count) + " groupName = " + *it + "\n"); count++; } }
            for (set<string>::iterator it = uniqueNames.begin(); it != uniqueNames.end(); it++) {  sraGroupNames.push_back(*it); }
        }

    
        return sraGroupNames;
    }
    catch(exception& e) {
        m->errorOut(e, "Oligos", "getSRAGroupNames");
        exit(1);
    }
}
//***************************************************************************************************************

int Oligos::readOligos(){
	try {
		ifstream inOligos;
        Utils util; util.openInputFile(oligosfile, inOligos);
		
		string type, oligo, roligo, group;
        bool pfUsesNone = false; bool prUsesNone = false; bool bfUsesNone = false; bool brUsesNone = false;
		
		while(!inOligos.eof()){
            
			inOligos >> type;
            
		 	if (m->getDebug()) { m->mothurOut("[DEBUG]: reading type - " + type + ".\n"); }
            
			if(type[0] == '#'){
				while (!inOligos.eof())	{	char c = inOligos.get();  if (c == 10 || c == 13){	break;	}	} // get rest of line if there's any crap there
				util.gobble(inOligos);
			}
			else{
				util.gobble(inOligos);
				//make type case insensitive
				for(int i=0;i<type.length();i++){	type[i] = toupper(type[i]);  }
				
				inOligos >> oligo;
                
                if (m->getDebug()) { m->mothurOut("[DEBUG]: reading - " + oligo + ".\n"); }
				
				for(int i=0;i<oligo.length();i++){
					oligo[i] = toupper(oligo[i]);
					if(oligo[i] == 'U')	{	oligo[i] = 'T';	}
				}
				
				if(type == "FORWARD"){
					group = "";
					
					// get rest of line in case there is a primer name
					while (!inOligos.eof())	{
						char c = inOligos.get();
						if (c == 10 || c == 13 || c == -1){	break;	}
						else if (c == 32 || c == 9){;} //space or tab
						else { 	group += c;  }
					}
					
					//check for repeat barcodes
					map<string, int>::iterator itPrime = primers.find(oligo);
					if (itPrime != primers.end()) { m->mothurOut("[WARNING]: primer " + oligo + " is in your oligos file already, disregarding."); m->mothurOutEndLine();  }
                    else {
                        if (m->getDebug()) {  if (group != "") { m->mothurOut("[DEBUG]: reading group " + group + ".\n"); }else{ m->mothurOut("[DEBUG]: no group for primer " + oligo + ".\n"); }  }
                    
                        primers[oligo]=indexPrimer; indexPrimer++;
                        primerNameVector.push_back(group);
                    }
				}
                else if (type == "PRIMER"){
                    util.gobble(inOligos);
					
                    inOligos >> roligo;
                    
                    for(int i=0;i<roligo.length();i++){
                        roligo[i] = toupper(roligo[i]);
                        if(roligo[i] == 'U')	{	roligo[i] = 'T';	}
                    }
                    
                    if (oligo == "NONE")        {  pfUsesNone = true; }
                    else if (roligo == "NONE")  { prUsesNone = true;  }
                   
                    if (roligo != "NONE") {
                        if (reversePairs) {  roligo = reverseOligo(roligo); }
                    }
                    group = "";
                    
					// get rest of line in case there is a primer name
					while (!inOligos.eof())	{
						char c = inOligos.get();
						if (c == 10 || c == 13 || c == -1){	break;	}
						else if (c == 32 || c == 9){;} //space or tab
						else { 	group += c;  }
					}
                    
                    oligosPair newPrimer(oligo, roligo);
                    
                    if (m->getDebug()) { m->mothurOut("[DEBUG]: primer pair " + newPrimer.forward + " " + newPrimer.reverse + ", and group = " + group + ".\n"); }
					
					//check for repeat barcodes
                    string tempPair = oligo+roligo;
                    if (uniquePrimers.count(tempPair) != 0) { m->mothurOut("primer pair " + newPrimer.forward + " " + newPrimer.reverse + " is in your oligos file already, disregarding."); m->mothurOutEndLine();  }
                    else { uniquePrimers.insert(tempPair);
					
                        if (m->getDebug()) {  if (group != "") { m->mothurOut("[DEBUG]: reading group " + group + ".\n"); }else{ m->mothurOut("[DEBUG]: no group for primer pair " + newPrimer.forward + " " + newPrimer.reverse + ".\n"); }  }
                
                        pairedPrimers[indexPairedPrimer]=newPrimer; indexPairedPrimer++;
                        primerNameVector.push_back(group);
                        hasPPrimers = true;
                    }
                }
				else if(type == "REVERSE"){
                    string oligoRC = reverseOligo(oligo);
					revPrimer.push_back(oligoRC);
				}
				else if(type == "BARCODE"){
					inOligos >> group;
                    
                    //barcode lines can look like   BARCODE   atgcatgc   groupName  - for 454 seqs
                    //or                            BARCODE   atgcatgc   atgcatgc    groupName  - for illumina data that has forward and reverse info
                    
                    string temp = "";
                    while (!inOligos.eof())	{
						char c = inOligos.get();
						if (c == 10 || c == 13 || c == -1){	break;	}
						else if (c == 32 || c == 9){;} //space or tab
						else { 	temp += c;  }
					}
					
                    //then this is illumina data with 4 columns
                    if (temp != "") {
                        hasPBarcodes = true;
                        string reverseBarcode = group; //reverseOligo(group); //reverse barcode
                        group = temp;
                        
                        for(int i=0;i<reverseBarcode.length();i++){
                            reverseBarcode[i] = toupper(reverseBarcode[i]);
                            if(reverseBarcode[i] == 'U')	{	reverseBarcode[i] = 'T';	}
                        }
                        
                        if (oligo == "NONE")                {  bfUsesNone = true; }
                        else if (reverseBarcode == "NONE")  { brUsesNone = true;  }
                        
                        if (reverseBarcode != "NONE") {
                            if (reversePairs) {   reverseBarcode = reverseOligo(reverseBarcode); }
                        }
                        oligosPair newPair(oligo, reverseBarcode);
                        
                        if (m->getDebug()) { m->mothurOut("[DEBUG]: barcode pair " + newPair.forward + " " + newPair.reverse + ", and group = " + group + ".\n"); }
                        
                        //check for repeat barcodes
                        string tempPair = oligo+reverseBarcode;
                        if (uniqueBarcodes.count(tempPair) != 0) { m->mothurOut("barcode pair " + newPair.forward + " " + newPair.reverse +  " is in your oligos file already, disregarding."); m->mothurOutEndLine();  }
                        else {
                            uniqueBarcodes.insert(tempPair);
                            pairedBarcodes[indexPairedBarcode]=newPair; indexPairedBarcode++;
                            barcodeNameVector.push_back(group);
                        }
                    }else {
                        //check for repeat barcodes
                        map<string, int>::iterator itBar = barcodes.find(oligo);
                        if (itBar != barcodes.end()) { m->mothurOut("[WARNING]: barcode " + oligo + " is in your oligos file already, disregarding."); m->mothurOutEndLine();  }
                        else {
                            barcodes[oligo]=indexBarcode; indexBarcode++;
                            barcodeNameVector.push_back(group);
                        }
                    }
				}else if(type == "LINKER"){
					linker.push_back(oligo);
				}else if(type == "SPACER"){
					spacer.push_back(oligo);
				}
				else{	m->mothurOut("[WARNING]: " + type + " is not recognized as a valid type. Choices are forward, reverse, and barcode. Ignoring " + oligo + "."); m->mothurOutEndLine(); }
			}
			util.gobble(inOligos);
		}
		inOligos.close();
		
        if ((linker.size() == 0) && (spacer.size() == 0) && (pairedBarcodes.size() == 0) && (barcodes.size() == 0) && (pairedPrimers.size() == 0) && (primers.size() == 0) && (revPrimer.size() == 0)) { m->mothurOut("[ERROR]: invalid oligos file, quitting.\n"); m->setControl_pressed(true); return 0; }
        
        if (hasPBarcodes || hasPPrimers) {
            pairedOligos = true;
            if ((primers.size() != 0) || (barcodes.size() != 0) || (linker.size() != 0) || (spacer.size() != 0) || (revPrimer.size() != 0)) { m->setControl_pressed(true);  m->mothurOut("[ERROR]: cannot mix paired primers and barcodes with non paired or linkers and spacers, quitting.\n");  return 0; }
            
            //check for "NONE" to make sure if none is used then all primers in that position are NONE
            //ex. Can't have: PRIMER NONE reversePrimer and PRIMER fowardPrimer reversePrimer in same file
            if (bfUsesNone) {
                bool allNONE = true;
                for(map<int, oligosPair>::iterator itBar = pairedBarcodes.begin();itBar != pairedBarcodes.end();itBar++){
                    if ((itBar->second).forward != "NONE") {
                        allNONE = false;
                        break;
                    }
                }
                if (!allNONE) {
                    m->setControl_pressed(true);  m->mothurOut("[ERROR]: cannot mix forwardBarcode=NONE and forwardBarcode=barcodeString in same file. Mothur assumes all sequences have forward barcodes or all do not, quitting."); m->mothurOutEndLine();  return 0;
                }
            }
            
            if (brUsesNone) {
                bool allNONE = true;
                for(map<int, oligosPair>::iterator itBar = pairedBarcodes.begin();itBar != pairedBarcodes.end();itBar++){
                    if ((itBar->second).reverse != "NONE") {
                        allNONE = false;
                        break;
                    }
                }
                if (!allNONE) {
                    m->setControl_pressed(true);  m->mothurOut("[ERROR]: cannot mix reverseBarcode=NONE and reverseBarcode=barcodeString in same file. Mothur assumes all sequences have reverse barcodes or all do not, quitting."); m->mothurOutEndLine();  return 0;
                }
            }
            
            if (pfUsesNone) {
                bool allNONE = true;
                for(map<int, oligosPair>::iterator itPrimer = pairedPrimers.begin();itPrimer != pairedPrimers.end(); itPrimer++){
                    if ((itPrimer->second).forward != "NONE") {
                        allNONE = false;
                        break;
                    }
                }
                if (!allNONE) {
                    m->setControl_pressed(true);  m->mothurOut("[ERROR]: cannot mix forwardPrimer=NONE and forwardPrimer=primerString in same file. Mothur assumes all sequences have forward primers or all do not, quitting."); m->mothurOutEndLine();  return 0;
                }
            }
            
            if (prUsesNone) {
                bool allNONE = true;
                for(map<int, oligosPair>::iterator itPrimer = pairedPrimers.begin();itPrimer != pairedPrimers.end(); itPrimer++){
                    if ((itPrimer->second).reverse != "NONE") {
                        allNONE = false;
                        break;
                    }
                }
                if (!allNONE) {
                    m->setControl_pressed(true);  m->mothurOut("[ERROR]: cannot mix reversePrimer=NONE and reversePrimer=primerString in same file. Mothur assumes all sequences have reverse primers or all do not, quitting."); m->mothurOutEndLine();  return 0;
                }
            }
        }
        
        
		//add in potential combos
		if(barcodeNameVector.size() == 0){
            if (pairedOligos) {
                oligosPair newPair("", "");
                pairedBarcodes[0] = newPair;
            }else {
                barcodes[""] = 0;
            }
			barcodeNameVector.push_back("");
		}
		
		if(primerNameVector.size() == 0){
            if (pairedOligos) {
                oligosPair newPair("", "");
                pairedPrimers[0] = newPair;
            }else {
                primers[""] = 0;
            }
			primerNameVector.push_back("");
		}
		
		
        if (pairedOligos) {
            for(map<int, oligosPair>::iterator itBar = pairedBarcodes.begin();itBar != pairedBarcodes.end();itBar++){
                for(map<int, oligosPair>::iterator itPrimer = pairedPrimers.begin();itPrimer != pairedPrimers.end(); itPrimer++){
                    
                    string primerName = primerNameVector[itPrimer->first];
                    string barcodeName = barcodeNameVector[itBar->first];
                    
                    if (m->getDebug()) {  m->mothurOut("[DEBUG]: primerName = " + primerName + " barcodeName = " + barcodeName + "\n");  }
                    
                    if ((primerName == "ignore") || (barcodeName == "ignore")) { if (m->getDebug()) {  m->mothurOut("[DEBUG]: in ignore. \n");  }  } //do nothing
                    else if ((primerName == "") && (barcodeName == "")) { if (m->getDebug()) {  m->mothurOut("[DEBUG]: in blank. \n");  }  } //do nothing
                    else {
                        string comboGroupName = "";
                        string fastqFileName = "";
                        
                        if(primerName == ""){
                            comboGroupName = barcodeNameVector[itBar->first];
                        }else{
                            if(barcodeName == ""){
                                comboGroupName = primerNameVector[itPrimer->first];
                            }
                            else{
                                comboGroupName = barcodeNameVector[itBar->first] + "." + primerNameVector[itPrimer->first];
                            }
                        }
                        
                        if (m->getDebug()) {  m->mothurOut("[DEBUG]: comboGroupName = " + comboGroupName +  "\n");  }
                        
                        uniqueNames.insert(comboGroupName);
                        
                        map<string, vector<string> >::iterator itGroup2Barcode = Group2Barcode.find(comboGroupName);
                        if (itGroup2Barcode == Group2Barcode.end()) {
                            vector<string> tempBarcodes; tempBarcodes.push_back((itBar->second).forward+"."+(itBar->second).reverse);
                            Group2Barcode[comboGroupName] = tempBarcodes;
                        }else {
                            Group2Barcode[comboGroupName].push_back((itBar->second).forward+"."+(itBar->second).reverse);
                        }
                        
                        itGroup2Barcode = Group2Primer.find(comboGroupName);
                        if (itGroup2Barcode == Group2Primer.end()) {
                            vector<string> tempPrimers; tempPrimers.push_back((itPrimer->second).forward+"."+(itPrimer->second).reverse);
                            Group2Primer[comboGroupName] = tempPrimers;
                        }else {
                            Group2Primer[comboGroupName].push_back((itPrimer->second).forward+"."+(itPrimer->second).reverse);
                        }
                    }
                }
            }
        }else {
            for(map<string, int>::iterator itBar = barcodes.begin();itBar != barcodes.end();itBar++){
                for(map<string, int>::iterator itPrimer = primers.begin();itPrimer != primers.end(); itPrimer++){
                    
                    string primerName = primerNameVector[itPrimer->second];
                    string barcodeName = barcodeNameVector[itBar->second];
                    
                    if ((primerName == "ignore") || (barcodeName == "ignore")) { } //do nothing
                    else if ((primerName == "") && (barcodeName == "")) { } //do nothing
                    else {
                        string comboGroupName = "";
                        string fastqFileName = "";
                        
                        if(primerName == ""){
                            comboGroupName = barcodeNameVector[itBar->second];
                        }
                        else{
                            if(barcodeName == ""){
                                comboGroupName = primerNameVector[itPrimer->second];
                            }
                            else{
                                comboGroupName = barcodeNameVector[itBar->second] + "." + primerNameVector[itPrimer->second];
                            }
                        }
                        uniqueNames.insert(comboGroupName);
                        
                        map<string, vector<string> >::iterator itGroup2Barcode = Group2Barcode.find(comboGroupName);
                        if (itGroup2Barcode == Group2Barcode.end()) {
                            vector<string> tempBarcodes; tempBarcodes.push_back(itBar->first);
                            Group2Barcode[comboGroupName] = tempBarcodes;
                        }else {
                            Group2Barcode[comboGroupName].push_back(itBar->first);
                        }
                        
                        itGroup2Barcode = Group2Primer.find(comboGroupName);
                        if (itGroup2Barcode == Group2Primer.end()) {
                            vector<string> tempPrimers; tempPrimers.push_back(itPrimer->first);
                            Group2Primer[comboGroupName] = tempPrimers;
                        }else {
                            Group2Primer[comboGroupName].push_back(itPrimer->first);
                        }
                    }
                }
            }
        }
        
        
        if (m->getDebug()) { int count = 0; for (set<string>::iterator it = uniqueNames.begin(); it != uniqueNames.end(); it++) { m->mothurOut("[DEBUG]: " + toString(count) + " groupName = " + *it + "\n"); count++; } }
        
        Groups.clear();
        for (set<string>::iterator it = uniqueNames.begin(); it != uniqueNames.end(); it++) {  Groups.push_back(*it); }
        
        return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "Oligos", "readOligos");
		exit(1);
	}
}
//********************************************************************/
vector<string> Oligos::getBarcodes(string groupName){
	try {
        vector<string> thisGroupsBarcodes;
        
        map<string, vector<string> >::iterator it = Group2Barcode.find(groupName);
        
        if (it == Group2Barcode.end()) {  m->mothurOut("[ERROR]: no barcodes found for group " + groupName + ".\n"); m->setControl_pressed(true);
        }else { thisGroupsBarcodes = it->second; }
        
        return thisGroupsBarcodes;
    }
	catch(exception& e) {
		m->errorOut(e, "Oligos", "getBarcodes");
		exit(1);
	}
}
//********************************************************************/
vector<string> Oligos::getPrimers(string groupName){
	try {
        vector<string> thisGroupsPrimers;
        
        map<string, vector<string> >::iterator it = Group2Primer.find(groupName);
        
        if (it == Group2Primer.end()) {  m->mothurOut("[ERROR]: no primers found for group " + groupName + ".\n"); m->setControl_pressed(true);
        }else { thisGroupsPrimers = it->second; }
        
        return thisGroupsPrimers;
    }
	catch(exception& e) {
		m->errorOut(e, "Oligos", "getPrimers");
		exit(1);
	}
}
//********************************************************************/
//can't have paired and unpaired so this function will either run the paired map or the unpaired
map<int, oligosPair> Oligos::getReorientedPairedPrimers(){
	try {
        map<int, oligosPair> rpairedPrimers;
        
        for (map<int, oligosPair>::iterator it = pairedPrimers.begin(); it != pairedPrimers.end(); it++) {
            string forward = (it->second).reverse;
            if (reversePairs) { forward = reverseOligo(forward); }
            string reverse = (it->second).forward;
            if (reversePairs) { reverse = reverseOligo(reverse); }
            oligosPair tempPair(forward, reverse); //reversePrimer, rc ForwardPrimer
            rpairedPrimers[it->first] = tempPair;
        }
        
        
        for (map<string, int>::iterator it = primers.begin(); it != primers.end(); it++) {
            oligosPair tempPair("", reverseOligo((it->first))); //reverseBarcode, rc ForwardBarcode
            rpairedPrimers[it->second] = tempPair;
        }
        
        return rpairedPrimers;
    }
	catch(exception& e) {
		m->errorOut(e, "Oligos", "getReorientedPairedPrimers");
		exit(1);
	}
}
//********************************************************************/
//can't have paired and unpaired so this function will either run the paired map or the unpaired
map<int, oligosPair> Oligos::getReorientedPairedBarcodes(){
	try {
        map<int, oligosPair> rpairedBarcodes;
        
        for (map<int, oligosPair>::iterator it = pairedBarcodes.begin(); it != pairedBarcodes.end(); it++) {
            string forward = (it->second).reverse;
            if (reversePairs) { forward = reverseOligo(forward); }
            string reverse = (it->second).forward;
            if (reversePairs) { reverse = reverseOligo(reverse); }
            oligosPair tempPair(forward, reverse); //reversePrimer, rc ForwardPrimer
            rpairedBarcodes[it->first] = tempPair;
        }
        
        for (map<string, int>::iterator it = barcodes.begin(); it != barcodes.end(); it++) {
            oligosPair tempPair("", reverseOligo((it->first))); //reverseBarcode, rc ForwardBarcode
            rpairedBarcodes[it->second] = tempPair;
        }

        return rpairedBarcodes;
    }
	catch(exception& e) {
		m->errorOut(e, "Oligos", "getReorientedPairedBarcodes");
		exit(1);
	}
}

//********************************************************************/
string Oligos::reverseOligo(string oligo){
	try {
        
        if (oligo == "NONE") { return "NONE"; }
        
        string reverse = "";
        
        for(int i=oligo.length()-1;i>=0;i--){
            
            if(oligo[i] == 'A')		{	reverse += 'T';	}
            else if(oligo[i] == 'T'){	reverse += 'A';	}
            else if(oligo[i] == 'U'){	reverse += 'A';	}
            
            else if(oligo[i] == 'G'){	reverse += 'C';	}
            else if(oligo[i] == 'C'){	reverse += 'G';	}
            
            else if(oligo[i] == 'R'){	reverse += 'Y';	}
            else if(oligo[i] == 'Y'){	reverse += 'R';	}
            
            else if(oligo[i] == 'M'){	reverse += 'K';	}
            else if(oligo[i] == 'K'){	reverse += 'M';	}
            
            else if(oligo[i] == 'W'){	reverse += 'W';	}
            else if(oligo[i] == 'S'){	reverse += 'S';	}
            
            else if(oligo[i] == 'B'){	reverse += 'V';	}
            else if(oligo[i] == 'V'){	reverse += 'B';	}
            
            else if(oligo[i] == 'D'){	reverse += 'H';	}
            else if(oligo[i] == 'H'){	reverse += 'D';	}
            
            else						{	reverse += 'N';	}
        }
        
        
        return reverse;
    }
	catch(exception& e) {
		m->errorOut(e, "Oligos", "reverseOligo");
		exit(1);
	}
}
//********************************************************************/
string Oligos::getBarcodeName(int index){
	try {
        string name = "";
        
        if ((index >= 0) && (index < barcodeNameVector.size())) { name = barcodeNameVector[index]; }
            
        return name;
    }
	catch(exception& e) {
		m->errorOut(e, "Oligos", "getBarcodeName");
		exit(1);
	}
}
//********************************************************************/
string Oligos::getPrimerName(int index){
    try {
        string name = "";
        
        if ((index >= 0) && (index < primerNameVector.size())) { name = primerNameVector[index]; }
        
        return name;
    }
    catch(exception& e) {
        m->errorOut(e, "Oligos", "getPrimerName");
        exit(1);
    }
}
//********************************************************************/
string Oligos::getGroupName(int barcodeIndex, int primerIndex){
    try {
        
        string thisGroup = "";
            if(numBarcodes != 0){
                thisGroup = getBarcodeName(barcodeIndex);
                if (numFPrimers != 0) {
                    if (getPrimerName(primerIndex) != "") {
                        if(thisGroup != "") {
                            thisGroup += "." + getPrimerName(primerIndex);
                        }else {
                            thisGroup = getPrimerName(primerIndex);
                        }
                    } 
                }
            }
        
        return thisGroup;
    }
    catch(exception& e) {
        m->errorOut(e, "Oligos", "getGroupName");
        exit(1);
    }
}

/**************************************************************************************************/

