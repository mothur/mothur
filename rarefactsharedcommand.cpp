/*
 *  rarefactsharedcommand.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/6/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "rarefactsharedcommand.h"
#include "sharedsobs.h"

//**********************************************************************************************************************

RareFactSharedCommand::RareFactSharedCommand(){
	try {
		globaldata = GlobalData::getInstance();
		string fileNameRoot;
		fileNameRoot = getRootName(globaldata->inputFileName);
		format = globaldata->getFormat();
		validCalculator = new ValidCalculators();
				
		int i;
		for (i=0; i<globaldata->Estimators.size(); i++) {
			if (validCalculator->isValidCalculator("sharedrarefaction", globaldata->Estimators[i]) == true) { 
				if (globaldata->Estimators[i] == "sharedobserved") { 
					rDisplays.push_back(new RareDisplay(new SharedSobs(), new SharedThreeColumnFile(fileNameRoot+"shared.rarefaction", "")));
				}
			}
		}
		
		//reset calc for next command
		globaldata->setCalc("");

	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the RareFactSharedCommand class Function RareFactSharedCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the RareFactSharedCommand class function RareFactSharedCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
			
}

//**********************************************************************************************************************

RareFactSharedCommand::~RareFactSharedCommand(){
	delete order;
	delete input;
	delete rCurve;
	delete read;
}

//**********************************************************************************************************************

int RareFactSharedCommand::execute(){
	try {
		int count = 1;
		
		//if the users entered no valid calculators don't execute command
		if (rDisplays.size() == 0) { return 0; }

		if (format == "sharedfile") {
			read = new ReadPhilFile(globaldata->inputFileName);	
			read->read(&*globaldata); 
			
			input = globaldata->ginput;
			order = input->getSharedOrderVector();
		}else {
			//you are using a list and a groupfile
			read = new ReadPhilFile(globaldata->inputFileName);	
			read->read(&*globaldata); 
		
			input = globaldata->ginput;
			SharedList = globaldata->gSharedList;
			order = SharedList->getSharedOrderVector();
		}
		
		//set users groups
		setGroups();
		
		while(order != NULL){
		
			if(globaldata->allLines == 1 || globaldata->lines.count(count) == 1 || globaldata->labels.count(order->getLabel()) == 1){
				//create collectors curve
				rCurve = new Rarefact(order, rDisplays);
				convert(globaldata->getFreq(), freq);
				convert(globaldata->getIters(), nIters);
				rCurve->getSharedCurve(freq, nIters);
			
				delete rCurve;
			
				cout << order->getLabel() << '\t' << count << endl;
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
	
		for(int i=0;i<rDisplays.size();i++){	delete rDisplays[i];	}	
		
		//reset groups parameter
		globaldata->Groups.clear();  globaldata->setGroups("");

		return 0;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the RareFactSharedCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the RareFactSharedCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}


//**********************************************************************************************************************

void RareFactSharedCommand::setGroups() {
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
				if ((globaldata->Groups.size() == 0) || (globaldata->Groups.size() == 1)) { 
					cout << "When using the groups parameter you must have at least 2 valid groups. I will run the command using all the groups in your groupfile." << endl; 
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
		cout << "Standard Error: " << e.what() << " has occurred in the RareFactSharedCommand class Function setGroups. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the RareFactSharedCommand class function setGroups. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		

}
/***********************************************************/

