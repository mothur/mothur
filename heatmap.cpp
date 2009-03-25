/*
 *  heatmap.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 3/25/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "heatmap.h"

//**********************************************************************************************************************
HeatMap::HeatMap(){
	try {
		globaldata = GlobalData::getInstance();
		format = globaldata->getFormat();
		
		if (format != "list") {	 setGroups();  }
		
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the HeatMap class Function HeatMap. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the HeatMap class function HeatMap. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
//**********************************************************************************************************************
void HeatMap::getPic(OrderVector* order) {
	try {
		sabund = order->getSAbundVector();
		string filename = getRootName(globaldata->inputFileName) + order->getLabel();
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the HeatMap class Function getPic. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the HeatMap class function getPic. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
//**********************************************************************************************************************
void HeatMap::getPic(SharedOrderVector* sharedorder) {
	try {
		getSharedVectors(sharedorder);
		
		
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the HeatMap class Function getPic. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the HeatMap class function getPic. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
//**********************************************************************************************************************
void HeatMap::getSharedVectors(SharedOrderVector* order){
	try {
		lookup.clear();
		
		vector<SharedRAbundVector*> templookup;
		
		//create and initialize vector of sharedvectors, one for each group
		for (int i = 0; i < globaldata->Groups.size(); i++) { 
			SharedRAbundVector* temp = new SharedRAbundVector(order->getNumBins());
			temp->setLabel(order->getLabel());
			temp->setGroup(globaldata->Groups[i]);
			templookup.push_back(temp);
		}
		
		int numSeqs = order->size();
		//sample all the members
		for(int i=0;i<numSeqs;i++){
			//get first sample
			individual chosen = order->get(i);
			int abundance; 
					
			//set info for sharedvector in chosens group
			for (int j = 0; j < templookup.size(); j++) { 
				if (chosen.group == templookup[j]->getGroup()) {
					 abundance = templookup[j]->getAbundance(chosen.bin);
					 templookup[j]->set(chosen.bin, (abundance + 1), chosen.group);
					 break;
				}
			}
		}
		
		//convert templookups rabunds to lookups sabunds
		for (int j = 0; j < templookup.size(); j++) { 
			lookup.push_back(templookup[j]->getSharedSAbundVector());
			delete templookup[j];
		}
		
		
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the HeatMap class Function getSharedVectors. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the HeatMap class function getSharedVectors. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}

}

//**********************************************************************************************************************
void HeatMap::setGroups() {
	try {
		//if the user has not entered specific groups to analyze then do them all
		if (globaldata->Groups.size() != 0) {
			if (globaldata->Groups[0] != "all") {
				//check that groups are valid
				for (int i = 0; i < globaldata->Groups.size(); i++) {
					if (globaldata->gGroupmap->isValidGroup(globaldata->Groups[i]) != true) {
						cout << globaldata->Groups[i] << " is not a valid group, and will be disregarded." << endl;
						// erase the invalid group from globaldata->Groups
						globaldata->Groups.erase(globaldata->Groups.begin()+i);
					}
				}
			
				//if the user only entered invalid groups
				if (globaldata->Groups.size() == 0) { 
					cout << "When using the groups parameter you must have at least 1 valid groups. I will run the command using all the groups in your groupfile." << endl; 
					for (int i = 0; i < globaldata->gGroupmap->namesOfGroups.size(); i++) {
						globaldata->Groups.push_back(globaldata->gGroupmap->namesOfGroups[i]);
					}
				}
			}else{//user has enter "all" and wants the default groups
				globaldata->Groups.clear();
				for (int i = 0; i < globaldata->gGroupmap->namesOfGroups.size(); i++) {
					globaldata->Groups.push_back(globaldata->gGroupmap->namesOfGroups[i]);
				}
				globaldata->setGroups("");
			}
		}else {
			for (int i = 0; i < globaldata->gGroupmap->namesOfGroups.size(); i++) {
				globaldata->Groups.push_back(globaldata->gGroupmap->namesOfGroups[i]);
			}
		}
		
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the HeatMap class Function setGroups. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the HeatMap class function setGroups. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		

}
/***********************************************************/
