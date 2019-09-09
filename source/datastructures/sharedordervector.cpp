/*
 *  sharedSharedOrderVector.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 12/9/08.
 *  Copyright 2008 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "sharedordervector.h"


/***********************************************************************/

SharedOrderVector::SharedOrderVector() : DataVector(), maxRank(0), numBins(0), numSeqs(0)  {}

/***********************************************************************/

SharedOrderVector::SharedOrderVector(string id, vector<individual>  ov) : 
											DataVector(id), data(ov)
{
		updateStats();
}

/***********************************************************************/
//This function is used to read a .shared file for the collect.shared, rarefaction.shared and summary.shared commands
//if you don't use a list and groupfile.  

SharedOrderVector::SharedOrderVector(ifstream& f, vector<string>& userGroups, string& previousLabel) : DataVector() {  //reads in a shared file
	try {
		maxRank = 0; numBins = 0; numSeqs = 0;
        int numUserGroups = userGroups.size();
				
		int num, inputData;
		numSeqs = 0;
		string holdLabel, nextLabel, groupN;
		individual newguy;
		
		//read in first row since you know there is at least 1 group.
		//are we at the beginning of the file??
		if (previousLabel == "") {
			f >> label; 
			
			//is this a shared file that has headers
			if (label == "label") { 
				//gets "group"
				f >> label; util.gobble(f);
				
				//gets "numOtus"
				f >> label; util.gobble(f);
				
				//eat rest of line
				label = util.getline(f); util.gobble(f);
				
				//parse labels to save
				istringstream iStringStream(label);
				while(!iStringStream.eof()){
					if (m->getControl_pressed()) { break; }
					string temp;
					iStringStream >> temp;  util.gobble(iStringStream);
					
					currentLabels.push_back(temp);
				}
				
				f >> label;
			}
		}else { label = previousLabel; }
		
		//read in first row since you know there is at least 1 group.
		f >> groupN >> num;
        bool readData = false;
        if (numUserGroups == 0) { //user has not specified groups, so we will use all of them
            userGroups.push_back(groupN);
            readData = true;
        }else{
            if (util.inUsersGroups(groupN, userGroups)) { readData = true; }
            //else - skipline because you are a group we dont care about
        }
		
		holdLabel = label;
        
        if (readData) {
            //save group in groupmap
            setNamesOfGroups(groupN);
            
            for(int i=0;i<num;i++){ //for each otu
                f >> inputData; //abundance of the otu
                
                for (int j = 0; j < inputData; j++) { //for each abundance
                    push_back(i, groupN);
                    numSeqs++;
                }
            }
        } else { util.getline(f); }
		
		util.gobble(f);
		
		if (!(f.eof())) { f >> nextLabel; }
		
		//read the rest of the groups info in
		while ((nextLabel == holdLabel) && (f.eof() != true)) {
			f >> groupN >> num;
            
            bool readData = false;
            if (numUserGroups == 0) { //user has not specified groups, so we will use all of them
                userGroups.push_back(groupN);
                readData = true;
            }else{
                if (util.inUsersGroups(groupN, userGroups)) { readData = true; }
                //else - skipline because you are a group we dont care about
            }

            if (readData) {
                //save group in groupmap
                setNamesOfGroups(groupN);
                
                for(int i=0;i<num;i++){
                    f >> inputData;
                    
                    for (int j = 0; j < inputData; j++) {
                        push_back(i, groupN);
                        numSeqs++;
                    }
                }
            }else { util.getline(f); }
			
			util.gobble(f);
				
			if (f.eof() != true) { f >> nextLabel; }

		}
		
		previousLabel = nextLabel;
        
        sort(userGroups.begin(), userGroups.end());
        for (int i = 0; i < userGroups.size(); i++) { setNamesOfGroups(userGroups[i]); }
        
		updateStats();
        
	}
	catch(exception& e) {
		m->errorOut(e, "SharedOrderVector", "SharedOrderVector");
		exit(1);
	}
}
/***********************************************************************/

int SharedOrderVector::getNumBins(){
	return numBins;
}

/***********************************************************************/

int SharedOrderVector::getNumSeqs(){
	return numSeqs;
}

/***********************************************************************/

int SharedOrderVector::getMaxRank(){
	return maxRank;
}


/***********************************************************************/
void SharedOrderVector::set(int index, int binNumber, int abund, string groupName){
    setNamesOfGroups(groupName);
	data[index].group = groupName;
	data[index].binNumber = binNumber;
	//if (abund > maxRank) { maxRank = abund; }
	updateStats();
}

/***********************************************************************/

individual SharedOrderVector::get(int index){
	return data[index];			
}

/************************************************************/
void SharedOrderVector::setNamesOfGroups(string seqGroup) {
    int i, count;
    count = 0;
    for (i=0; i<allGroups.size(); i++) {
        if (allGroups[i] != seqGroup) {
            count++; //you have not found this group
        }else {
            break; //you already have it
        }
    }
    if (count == allGroups.size()) {
        allGroups.push_back(seqGroup); //new group
    }
}
/***********************************************************************/
//commented updateStats out to improve speed, but whoever calls this must remember to update when they are done with all the pushbacks they are doing 
void SharedOrderVector::push_back(int bin, string groupName){
	setNamesOfGroups(groupName);
    individual newGuy;
	newGuy.group = groupName;
	newGuy.binNumber = bin;
	data.push_back(newGuy);
}

/***********************************************************************/

void SharedOrderVector::print(ostream& output){
	try {
		output << label << '\t' << numSeqs;
	
		for(int i=0;i<data.size();i++){
			output << '\t' << data[i].binNumber;
		}
		output << endl;
	}
	catch(exception& e) {
		m->errorOut(e, "SharedOrderVector", "print");
		exit(1);
	}
}

/***********************************************************************/

void SharedOrderVector::clear(){
	numBins = 0;
	maxRank = 0;
	numSeqs = 0;
	data.clear();
}
/***********************************************************************/

void SharedOrderVector::resize(int){
	m->mothurOut("resize() did nothing in class SharedOrderVector");
}

/***********************************************************************/


vector<individual>::iterator SharedOrderVector::begin(){
	return data.begin();	
}

/***********************************************************************/

vector<individual>::iterator SharedOrderVector::end(){
	return data.end();		
}

/***********************************************************************/

int SharedOrderVector::size(){
	return data.size();					
}

/***********************************************************************/

RAbundVector SharedOrderVector::getRAbundVector(){
	try {
		RAbundVector rav(data.size());
	
		for(int i=0;i<numSeqs;i++){
			rav.set(data[i].binNumber, rav.get(data[i].binNumber) + 1);
		}	
		sort(rav.rbegin(), rav.rend());
		for(int i=numSeqs-1;i>=0;i--){
			if(rav.get(i) == 0){	rav.pop_back();	}
			else{ break; }
		}
		rav.setLabel(label);

		return rav;
	}
	catch(exception& e) {
		m->errorOut(e, "SharedOrderVector", "getRAbundVector");
		exit(1);
	}
}
/***********************************************************************/

OrderVector SharedOrderVector::getOrderVector(map<string,int>* nameMap = NULL) {
	try {
		OrderVector ov;
	
		for (int i = 0; i < data.size(); i++) {
			ov.push_back(data[i].binNumber);
		}
		
		util.mothurRandomShuffle(ov);

		ov.setLabel(label);	
		return ov;
	}
	catch(exception& e) {
		m->errorOut(e, "SharedOrderVector", "getOrderVector");
		exit(1);
	}
}


/***********************************************************************/

SAbundVector SharedOrderVector::getSAbundVector(){
	
	RAbundVector rav(this->getRAbundVector());
	return rav.getSAbundVector();

}
/***********************************************************************/
SharedRAbundVectors* SharedOrderVector::getSharedRAbundVector(string group) {
	try {
		SharedRAbundVector* sharedRav = new SharedRAbundVector(data.size());
		
		sharedRav->setLabel(label);
		sharedRav->setGroup(group);
		
		for (int i = 0; i < data.size(); i++) {
			if (data[i].group == group) {
				sharedRav->set(data[i].binNumber, sharedRav->get(data[i].binNumber) + 1);
			}
		}
        
        SharedRAbundVectors* lookup = new SharedRAbundVectors();
        lookup->setOTUNames(currentLabels);
        lookup->push_back(sharedRav);
		return lookup;
	}
	catch(exception& e) {
		m->errorOut(e, "SharedOrderVector", "getSharedRAbundVector");
		exit(1);
	}
}
/***********************************************************************/
SharedRAbundVectors* SharedOrderVector::getSharedRAbundVector() {
	try {
        SharedRAbundVectors* lookup = new SharedRAbundVectors();
        sort(allGroups.begin(), allGroups.end());
        
        //create and initialize vector of sharedvectors, one for each group
        for (int i = 0; i < allGroups.size(); i++) {
            SharedRAbundVector* temp = new SharedRAbundVector(numBins);
            temp->setLabel(getLabel());
            temp->setGroup(allGroups[i]);
            lookup->push_back(temp);
        }
        
        int numSeqs = size();
        //sample all the members
        for(int i=0;i<numSeqs;i++){
            //get first sample
            individual chosen = get(i);
            int abundance = lookup->get(chosen.binNumber, chosen.group);
            lookup->set(chosen.binNumber, (abundance + 1), chosen.group);
        }
        
        lookup->setOTUNames(currentLabels);
		
		return lookup;
	}
	catch(exception& e) {
		m->errorOut(e, "SharedOrderVector", "getSharedRAbundVector");
		exit(1);
	}
}
/***********************************************************************/

SharedOrderVector SharedOrderVector::getSharedOrderVector(){
	util.mothurRandomShuffle(*this);
	return *this;			
}

/***********************************************************************/

void SharedOrderVector::updateStats(){
	try {
		needToUpdate = 0;
		numSeqs = 0;
		numBins = 0;
		maxRank = 0;
	
		numSeqs = data.size();
        
        std::set<int> uniqueBins;
		for(int i=0;i<numSeqs;i++){
            uniqueBins.insert(data[i].binNumber);
			if(data[i].binNumber > maxRank)		{	maxRank = data[i].binNumber;		}
		}
        numBins = uniqueBins.size();
	}
	catch(exception& e) {
		m->errorOut(e, "SharedOrderVector", "updateStats");
		exit(1);
	}
}

/***********************************************************************/


