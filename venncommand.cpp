/*
 *  venncommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 3/30/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "venncommand.h"
#include "ace.h"
#include "sobs.h"
#include "chao1.h"
#include "jackknife.h"
#include "sharedsobscollectsummary.h"
#include "sharedchao1.h"
#include "sharedace.h"


//**********************************************************************************************************************

VennCommand::VennCommand(){
	try {
		globaldata = GlobalData::getInstance();
		format = globaldata->getFormat();
		validCalculator = new ValidCalculators();
		
		int i;
		
		if (format == "list") {
			for (i=0; i<globaldata->Estimators.size(); i++) {
				if (validCalculator->isValidCalculator("vennsingle", globaldata->Estimators[i]) == true) { 
					if (globaldata->Estimators[i] == "sobs") { 
						vennCalculators.push_back(new Sobs());
					}else if (globaldata->Estimators[i] == "chao") { 
						vennCalculators.push_back(new Chao1());
					}else if (globaldata->Estimators[i] == "ace") {
						convert(globaldata->getAbund(), abund);
						if(abund < 5)
							abund = 10;
						vennCalculators.push_back(new Ace(abund));
					}else if (globaldata->Estimators[i] == "jack") { 	
						vennCalculators.push_back(new Jackknife());
					}
				}
			}
		}else {
			for (i=0; i<globaldata->Estimators.size(); i++) {
				if (validCalculator->isValidCalculator("vennshared", globaldata->Estimators[i]) == true) { 
					if (globaldata->Estimators[i] == "sharedsobs") { 
						vennCalculators.push_back(new SharedSobsCS());
					}else if (globaldata->Estimators[i] == "sharedchao") { 
						vennCalculators.push_back(new SharedChao1());
					}else if (globaldata->Estimators[i] == "sharedace") { 
						vennCalculators.push_back(new SharedAce());
					}
				}
			}
		}
		
		venn = new Venn();
		
		//reset calc for next command
		globaldata->setCalc("");

		
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the VennCommand class Function VennCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the VennCommand class function VennCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}
//**********************************************************************************************************************

VennCommand::~VennCommand(){
	delete input;
	delete read;
	delete venn;
}

//**********************************************************************************************************************

int VennCommand::execute(){
	try {
		int count = 1;	
		
		//if the users entered no valid calculators don't execute command
		if (vennCalculators.size() == 0) { return 0; }
		
		if (format == "sharedfile") {
			//you have groups
			read = new ReadPhilFile(globaldata->inputFileName);	
			read->read(&*globaldata); 
			
			input = globaldata->ginput;
			order = input->getSharedOrderVector();
		}else if (format == "shared") {
			//you are using a list and a groupfile
			read = new ReadPhilFile(globaldata->inputFileName);	
			read->read(&*globaldata); 
		
			input = globaldata->ginput;
			SharedList = globaldata->gSharedList;
			order = SharedList->getSharedOrderVector();
		}else if (format == "list") {
			//you are using just a list file and have only one group
			read = new ReadPhilFile(globaldata->inputFileName);	
			read->read(&*globaldata); 
		
			ordersingle = globaldata->gorder;
			input = globaldata->ginput;
		}

		
		if (format != "list") {	
			
			setGroups();
			
			while(order != NULL){
		
				if(globaldata->allLines == 1 || globaldata->lines.count(count) == 1 || globaldata->labels.count(order->getLabel()) == 1){			
	
					cout << order->getLabel() << '\t' << count << endl;
					venn->getPic(order, vennCalculators);

				}
						
				//get next line to process
				if (format == "sharedfile") {
					order = input->getSharedOrderVector();
				}else {
					//you are using a list and a groupfile
					SharedList = input->getSharedListVector(); //get new list vector to process
					if (SharedList != NULL) {
						order = SharedList->getSharedOrderVector(); //gets new order vector with group info.
					}else {
						break;
					}
				}
				count++;
			}
			
			//reset groups parameter
			globaldata->Groups.clear();  globaldata->setGroups("");
			
		}else{
			while(ordersingle != NULL){
		
				if(globaldata->allLines == 1 || globaldata->lines.count(count) == 1 || globaldata->labels.count(ordersingle->getLabel()) == 1){			
	
					cout << ordersingle->getLabel() << '\t' << count << endl;
					venn->getPic(ordersingle, vennCalculators);
					
				}
				
				ordersingle = (input->getOrderVector());
				count++;
			}
		}
		
		return 0;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the VennCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the VennCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		
}

//**********************************************************************************************************************
void VennCommand::setGroups() {
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
					if (globaldata->gGroupmap->namesOfGroups.size() > 4) {
						cout << "When using the groups parameter you must have at least 1 valid group. I will run the command using the first four groups in your groupfile." << endl; 
						for (int i = 0; i < 4; i++) {
							globaldata->Groups.push_back(globaldata->gGroupmap->namesOfGroups[i]);
						}
					}else {
						cout << "When using the groups parameter you must have at least 1 valid group. I will run the command using all the groups in your groupfile." << endl; 
						for (int i = 0; i < globaldata->gGroupmap->namesOfGroups.size(); i++) {
							globaldata->Groups.push_back(globaldata->gGroupmap->namesOfGroups[i]);
						}
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
		
		
		//check to make sure their are only 3 groups
		if (globaldata->Groups.size() > 4) {
			cout << "You may only use up to 4 groups at a time with this command.  I will choose the first four and disregard the rest." << endl;
			for (int i = 4; i < globaldata->Groups.size(); i++) {
				globaldata->Groups.erase(globaldata->Groups.begin()+i);
			}
		}
		
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the VennCommand class Function setGroups. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the VennCommand class function setGroups. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		

}
/***********************************************************/

