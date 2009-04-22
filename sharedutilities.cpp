/*
 *  sharedutilities.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 4/9/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "sharedutilities.h"
#include "sharedrabundvector.h"
#include "sharedordervector.h"

/**************************************************************************************************/

void SharedUtil::getSharedVectors(vector<string> Groups, vector<SharedRAbundVector*>& lookup, SharedOrderVector* order) {
	try {
	
		//delete each sharedrabundvector in lookup
		for (int j = 0; j < lookup.size(); j++) {
			delete lookup[j];
		}
		
		lookup.clear();
		
		//create and initialize vector of sharedvectors, one for each group
		for (int i = 0; i < Groups.size(); i++) { 
			SharedRAbundVector* temp = new SharedRAbundVector(order->getNumBins());
			temp->setLabel(order->getLabel());
			temp->setGroup(Groups[i]);
			lookup.push_back(temp);
		}
	
		int numSeqs = order->size();
		//sample all the members
		for(int i=0;i<numSeqs;i++){
			//get first sample
			individual chosen = order->get(i);
			int abundance; 
					
			//set info for sharedvector in chosens group
			for (int j = 0; j < lookup.size(); j++) { 
				if (chosen.group == lookup[j]->getGroup()) {
					 abundance = lookup[j]->getAbundance(chosen.bin);
					 lookup[j]->set(chosen.bin, (abundance + 1), chosen.group);
					 break;
				}
			}
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SharedUtil class Function getSharedVectors. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SharedUtil class function getSharedVectors. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
/**************************************************************************************************/

void SharedUtil::getSharedVectorswithReplacement(vector<string> Groups, vector<SharedRAbundVector*>& lookup, SharedOrderVector* order) {
	try {
	
		//delete each sharedrabundvector in lookup
		for (int j = 0; j < lookup.size(); j++) {
			delete lookup[j];
		}
		
		lookup.clear();
		
		//create and initialize vector of sharedvectors, one for each group
		for (int i = 0; i < Groups.size(); i++) { 
			SharedRAbundVector* temp = new SharedRAbundVector(order->getNumBins());
			temp->setLabel(order->getLabel());
			temp->setGroup(Groups[i]);
			lookup.push_back(temp);
		}
	
		int numSeqs = order->size();
		
		//sample all the members
		for(int i=0;i<numSeqs;i++){
			//get random number
			int random = int((float)(i+1) * (float)(rand()) / ((float)RAND_MAX+1.0));
			individual chosen = order->get(random);

			int abundance; 
			//set info for sharedvector in chosens group
			for (int j = 0; j < lookup.size(); j++) { 
				if (chosen.group == lookup[j]->getGroup()) {
					 abundance = lookup[j]->getAbundance(chosen.bin);
					 lookup[j]->set(chosen.bin, (abundance + 1), chosen.group);
					 break;
				}
			}
		}
		
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SharedUtil class Function getSharedVectorswithReplacement. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SharedUtil class function getSharedVectorswithReplacement. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/**************************************************************************************************/
//need to have mode because different commands require different number of valid groups
void SharedUtil::setGroups(vector<string>& userGroups, vector<string>& allGroups) {
	try {
		if (userGroups.size() != 0) {
			if (userGroups[0] != "all") {
				//check that groups are valid
				for (int i = 0; i < userGroups.size(); i++) {
					if (isValidGroup(userGroups[i], allGroups) != true) {
						cout << userGroups[i] << " is not a valid group, and will be disregarded." << endl;
						// erase the invalid group from userGroups
						userGroups.erase(userGroups.begin()+i);
						i--;
					}
				}
				
				//if the user only entered invalid groups
				if (userGroups.size() == 0) { 
					cout << "When using the groups parameter you must have at least 1 valid groups. I will run the command using all the groups in your groupfile." << endl; 
					for (int i = 0; i < allGroups.size(); i++) {
						userGroups.push_back(allGroups[i]);
					}
				}

			}else{//user has enter "all" and wants the default groups
				userGroups.clear();
				for (int i = 0; i < allGroups.size(); i++) {
					userGroups.push_back(allGroups[i]);
				}
			}
		}else { //the user has not entered groups
			for (int i = 0; i < allGroups.size(); i++) {
				userGroups.push_back(allGroups[i]);
			}
		}
			
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SharedUtil class Function setGroups. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SharedUtil class function setGroups. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}

}
/**************************************************************************************************/
//need to have mode because different commands require different number of valid groups
void SharedUtil::setGroups(vector<string>& userGroups, vector<string>& allGroups, string mode) {
	try {
		if (userGroups.size() != 0) {
			if (userGroups[0] != "all") {
				//check that groups are valid
				for (int i = 0; i < userGroups.size(); i++) {
					if (isValidGroup(userGroups[i], allGroups) != true) {
						cout << userGroups[i] << " is not a valid group, and will be disregarded." << endl;
						// erase the invalid group from userGroups
						userGroups.erase(userGroups.begin()+i);
						i--;
					}
				}

			}else{//user has enter "all" and wants the default groups
				userGroups.clear();
				for (int i = 0; i < allGroups.size(); i++) {
					userGroups.push_back(allGroups[i]);
				}
			}
		}else { //the user has not entered groups
			for (int i = 0; i < allGroups.size(); i++) {
				userGroups.push_back(allGroups[i]);
			}
		}
			
		if ((mode == "collect") || (mode == "rarefact") || (mode == "summary") || (mode == "treegroup")) {
				//if the user only entered invalid groups
				if ((userGroups.size() == 0) || (userGroups.size() == 1)) { 
					cout << "When using the groups parameter you must have at least 2 valid groups. I will run the command using all the groups in your groupfile." << endl; 
					for (int i = 0; i < allGroups.size(); i++) {
						userGroups.push_back(allGroups[i]);
					}
				}
		}
	
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SharedUtil class Function setGroups. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SharedUtil class function setGroups. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}

}


/**************************************************************************************/
//for parsimony and unifrac commands you set pairwise groups as well as an allgroups in calc
void SharedUtil::setGroups(vector<string>& userGroups, vector<string>& allGroups, string& label, int& numGroups, string mode){  //globaldata->Groups, your tree or group map, allgroups, mode
	try {
		numGroups = 0;
		label = "";
		//if the user has not entered specific groups to analyze then do them all
		if (userGroups.size() != 0) {
			if (userGroups[0] != "all") {
				//check that groups are valid
				for (int i = 0; i < userGroups.size(); i++) {
					if (isValidGroup(userGroups[i], allGroups) != true) {
						cout << userGroups[i] << " is not a valid group, and will be disregarded." << endl;
						// erase the invalid group from globaldata->Groups
						userGroups.erase(userGroups.begin()+i);
						i--;
					}
				}
			}else { //users wants all groups
				userGroups.clear();
				for (int i=0; i < allGroups.size(); i++) { 
					if (allGroups[i] != "xxx") {
						userGroups.push_back(allGroups[i]);
					}
				}
			}
		}else { //the user has not entered groups
			for (int i=0; i < allGroups.size(); i++) { 
				if (allGroups[i] != "xxx") {
					if (mode == "weighted") {
						userGroups.push_back(allGroups[i]);
					}else {
						numGroups = 1;
						label += allGroups[i] + "-";
					}
				}
			}
			//rip extra - off allgroups 
			label = label.substr(0, label.length()-1);
		}
		
		if (mode == "weighted") {
			//if the user only entered invalid groups
			if (userGroups.size() == 0) { 
				for (int i=0; i < allGroups.size(); i++) { 
					if (allGroups[i] != "xxx") {
						userGroups.push_back(allGroups[i]);
					}
				}
				cout << "When using the groups parameter you must have at least 2 valid groups. I will run the command using all the groups in your groupfile." << endl; 
			}else if (userGroups.size() == 1) { 
				cout << "When using the groups parameter you must have at least 2 valid groups. I will run the command using all the groups in your groupfile." << endl;
				userGroups.clear();
				for (int i=0; i < allGroups.size(); i++) { 
					if (allGroups[i] != "xxx") {
						userGroups.push_back(allGroups[i]);
					}
				}
			}
			numGroups = userGroups.size();
			
		}else if ((mode == "unweighted") || (mode == "parsimony")) {
				//if the user only entered invalid groups
				if ((userGroups.size() == 0) && (numGroups == 0)) { 
					cout << "When using the groups parameter you must have at least 1 valid group. I will run the command using all the groups in your groupfile." << endl; 
					for (int i = 0; i < allGroups.size(); i++) {
						if (allGroups[i] != "xxx") {
							userGroups.push_back(allGroups[i]);
						}
					}
				}
				
				if (numGroups != 1) { numGroups = userGroups.size(); }
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SharedUtil class Function setGroups. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SharedUtil class function setGroups. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
/**************************************************************************************/
void SharedUtil::getCombos(vector<string>& groupComb, vector<string> userGroups, int& numComp) { //groupcomb, globaldata->Groups, numcomb
	try {
		//calculate number of comparisons i.e. with groups A,B,C = AB, AC, BC = 3;
		numComp = 0;
		for (int i=0; i< userGroups.size(); i++) { 
			numComp += i; 
			for (int l = i+1; l < userGroups.size(); l++) {
				//set group comparison labels
				groupComb.push_back(userGroups[i] + "-" + userGroups[l]);
			}
		} 
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SharedUtil class Function getCombos. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SharedUtil class function getCombos. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
/**************************************************************************************/
bool SharedUtil::isValidGroup(string groupname, vector<string> groups) {
	try {
		for (int i = 0; i < groups.size(); i++) {
			if (groupname == groups[i]) { return true; }
		}
		
		return false;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SharedUtil class Function isValidGroup. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SharedUtil class function isValidGroup. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/**************************************************************************************/
void SharedUtil::updateGroupIndex(vector<string>& userGroups, map<string, int>& index) {
	try {
		index.clear();
		for (int i = 0; i < userGroups.size(); i++) {
			index[userGroups[i]] = i;
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SharedUtil class Function updateGroupIndex. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SharedUtil class function updateGroupIndex. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
