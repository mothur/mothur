/*
 *  sharedSharedOrderVector.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 12/9/08.
 *  Copyright 2008 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "sharedordervector.h"
#include "sharedutilities.h"

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

SharedOrderVector::SharedOrderVector(ifstream& f) : DataVector() {  //reads in a shared file
	try {
		maxRank = 0; numBins = 0; numSeqs = 0;
				
		groupmap = new GroupMap(); 
		
		int num, inputData, count;
		count = 0;  numSeqs = 0;
		string holdLabel, nextLabel, groupN;
		individual newguy;
		
		//read in first row since you know there is at least 1 group.
		//are we at the beginning of the file??
		if (m->saveNextLabel == "") {  
			f >> label; 
			
			//is this a shared file that has headers
			if (label == "label") { 
				//gets "group"
				f >> label; m->gobble(f);
				
				//gets "numOtus"
				f >> label; m->gobble(f);
				
				//eat rest of line
				label = m->getline(f); m->gobble(f);
				
				//parse labels to save
				istringstream iStringStream(label);
				m->binLabelsInFile.clear();
				while(!iStringStream.eof()){
					if (m->control_pressed) { break; }
					string temp;
					iStringStream >> temp;  m->gobble(iStringStream);
					
					m->binLabelsInFile.push_back(temp);
				}
				
				f >> label;
			}
		}else { label = m->saveNextLabel; }
		
		//reset labels, currentLabels may have gotten changed as otus were eliminated because of group choices or sampling
		m->currentBinLabels = m->binLabelsInFile;
		
		//read in first row since you know there is at least 1 group.
		f >> groupN >> num;
		
		holdLabel = label;
		
		
		vector<string> allGroups;
		//save group in groupmap
		allGroups.push_back(groupN);
		groupmap->groupIndex[groupN] = 0;
		
		
		for(int i=0;i<num;i++){
			f >> inputData;
			
			for (int j = 0; j < inputData; j++) {
				push_back(i, i, groupN);
				numSeqs++;
			}
		}
		
		m->gobble(f); 
		
		if (!(f.eof())) { f >> nextLabel; }
		
		//read the rest of the groups info in
		while ((nextLabel == holdLabel) && (f.eof() != true)) {
			f >> groupN >> num;
			count++;
			
			
			//save group in groupmap
			allGroups.push_back(groupN);
			groupmap->groupIndex[groupN] = count;
			
			
			for(int i=0;i<num;i++){
				f >> inputData;
				
				for (int j = 0; j < inputData; j++) {
					push_back(i, i, groupN);
					numSeqs++;
				}
			}
			
			m->gobble(f);
				
			if (f.eof() != true) { f >> nextLabel; }

		}
		
		m->saveNextLabel = nextLabel;
			
		groupmap->setNamesOfGroups(allGroups);
		m->setAllGroups(allGroups);
		
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
	
	data[index].group = groupName;
	data[index].bin = binNumber;
	data[index].abundance = abund;
	//if (abund > maxRank) { maxRank = abund; }
	updateStats();
}

/***********************************************************************/

individual SharedOrderVector::get(int index){
	return data[index];			
}


/***********************************************************************/
//commented updateStats out to improve speed, but whoever calls this must remember to update when they are done with all the pushbacks they are doing 
void SharedOrderVector::push_back(int binNumber, int abund, string groupName){
	individual newGuy;
	newGuy.group = groupName;
	newGuy.abundance = abund;
	newGuy.bin = binNumber;
	data.push_back(newGuy);
	//numSeqs++;
	//numBins++;
	//if (abund > maxRank) { maxRank = abund; }
	
	//updateStats();
}

/***********************************************************************/

void SharedOrderVector::print(ostream& output){
	try {
		output << label << '\t' << numSeqs << '\t';
	
		for(int i=0;i<data.size();i++){
			output << data[i].bin << '\t';
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
			rav.set(data[i].bin, rav.get(data[i].bin) + 1);
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
		m->errorOut(e, "SharedOrderVector", "getRAbundVector");
		exit(1);
	}
}
/***********************************************************************/

OrderVector SharedOrderVector::getOrderVector(map<string,int>* nameMap = NULL) {
	try {
		OrderVector ov;
	
		for (int i = 0; i < data.size(); i++) {
			ov.push_back(data[i].bin);
		}
		
		random_shuffle(ov.begin(), ov.end());

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
SharedRAbundVector SharedOrderVector::getSharedRAbundVector(string group) {
	try {
		SharedRAbundVector sharedRav(data.size());
		
		sharedRav.setLabel(label);
		sharedRav.setGroup(group);
		
		for (int i = 0; i < data.size(); i++) {
			if (data[i].group == group) {
				sharedRav.set(data[i].abundance, sharedRav.getAbundance(data[i].abundance) + 1, data[i].group);
			}
		}
		return sharedRav;
	}
	catch(exception& e) {
		m->errorOut(e, "SharedOrderVector", "getSharedRAbundVector");
		exit(1);
	}
}
/***********************************************************************/
vector<SharedRAbundVector*> SharedOrderVector::getSharedRAbundVector() {
	try {
		SharedUtil* util;
		util = new SharedUtil();
		vector<SharedRAbundVector*> lookup;
		
		vector<string> Groups = m->getGroups();
		vector<string> allGroups = m->getAllGroups();
		util->setGroups(Groups, allGroups);
		util->getSharedVectors(Groups, lookup, this);
		m->setGroups(Groups);
		m->setAllGroups(allGroups);
		
		return lookup;
	}
	catch(exception& e) {
		m->errorOut(e, "SharedOrderVector", "getSharedRAbundVector");
		exit(1);
	}
}
/***********************************************************************/
SharedSAbundVector SharedOrderVector::getSharedSAbundVector(string group) {
	try {
		
		SharedRAbundVector sharedRav(this->getSharedRAbundVector(group));
		return sharedRav.getSharedSAbundVector();
				
	}
	catch(exception& e) {
		m->errorOut(e, "SharedOrderVector", "getSharedSAbundVector");
		exit(1);
	}
}

/***********************************************************************/

SharedOrderVector SharedOrderVector::getSharedOrderVector(){
	random_shuffle(data.begin(), data.end());
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
				
		vector<int> hold(numSeqs, 0);
		for(int i=0;i<numSeqs;i++){
			hold[data[i].bin] = hold[data[i].bin]+1;
		}	
		
		for(int i=0;i<numSeqs;i++){
			if(hold[i] > 0)				{	numBins++;				}
			if(hold[i] > maxRank)		{	maxRank = hold[i];		}
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "SharedOrderVector", "updateStats");
		exit(1);
	}
}

/***********************************************************************/


