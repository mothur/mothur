/*
 *  collectsharedcommand.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "collectsharedcommand.h"
#include "sharedchao1.h"
#include "sharedace.h"
#include "sharedjabund.h"
#include "sharedsorabund.h"
#include "sharedjclass.h"
#include "sharedsorclass.h"
#include "sharedjest.h"
#include "sharedsorest.h"
#include "sharedthetayc.h"
#include "sharedthetan.h"


//**********************************************************************************************************************

CollectSharedCommand::CollectSharedCommand(){
	try {
		globaldata = GlobalData::getInstance();
		string fileNameRoot;
		fileNameRoot = getRootName(globaldata->inputFileName);
		groupmap = globaldata->gGroupmap;
		
		int i;
		for (i=0; i<globaldata->sharedEstimators.size(); i++) {
			if (globaldata->sharedEstimators[i] == "sharedChao") { 
				cDisplays.push_back(new CollectDisplay(new SharedChao1(), new SharedOneColumnFile(fileNameRoot+"sharedChao", groupmap->namesOfGroups)));
			}else if (globaldata->sharedEstimators[i] == "sharedAce") { 
				cDisplays.push_back(new CollectDisplay(new SharedAce(), new SharedOneColumnFile(fileNameRoot+"sharedAce", groupmap->namesOfGroups)));
			}else if (globaldata->sharedEstimators[i] == "sharedJabund") { 	
				cDisplays.push_back(new CollectDisplay(new SharedJAbund(), new SharedOneColumnFile(fileNameRoot+"SharedJabund", groupmap->namesOfGroups)));
			}else if (globaldata->sharedEstimators[i] == "sharedSorensonAbund") { 
				cDisplays.push_back(new CollectDisplay(new SharedSorAbund(), new SharedOneColumnFile(fileNameRoot+"SharedSorensonAbund", groupmap->namesOfGroups)));
			}else if (globaldata->sharedEstimators[i] == "sharedJclass") { 
				cDisplays.push_back(new CollectDisplay(new SharedJclass(), new SharedOneColumnFile(fileNameRoot+"SharedJclass", groupmap->namesOfGroups)));
			}else if (globaldata->sharedEstimators[i] == "sharedSorClass") { 
				cDisplays.push_back(new CollectDisplay(new SharedSorClass(), new SharedOneColumnFile(fileNameRoot+"SharedSorClass", groupmap->namesOfGroups)));
			}else if (globaldata->sharedEstimators[i] == "sharedJest") { 
				cDisplays.push_back(new CollectDisplay(new SharedJest(), new SharedOneColumnFile(fileNameRoot+"SharedJest", groupmap->namesOfGroups)));
			}else if (globaldata->sharedEstimators[i] == "sharedSorEst") { 
				cDisplays.push_back(new CollectDisplay(new SharedSorEst(), new SharedOneColumnFile(fileNameRoot+"SharedSorEst", groupmap->namesOfGroups)));
			}else if (globaldata->sharedEstimators[i] == "SharedThetaYC") { 
				cDisplays.push_back(new CollectDisplay(new SharedThetaYC(), new SharedOneColumnFile(fileNameRoot+"SharedThetaYC", groupmap->namesOfGroups)));
			}else if (globaldata->sharedEstimators[i] == "SharedThetaN") { 
				cDisplays.push_back(new CollectDisplay(new SharedThetaN(), new SharedOneColumnFile(fileNameRoot+"SharedThetaN", groupmap->namesOfGroups)));
			}
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the CollectSharedCommand class Function CollectSharedCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the CollectSharedCommand class function CollectSharedCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
			
}

//**********************************************************************************************************************

CollectSharedCommand::~CollectSharedCommand(){
	delete order;
	delete input;
	delete cCurve;
	delete read;
}

//**********************************************************************************************************************

int CollectSharedCommand::execute(){
	try {
		int count = 1;
		read = new ReadPhilFile(globaldata->inputFileName);	
		read->read(&*globaldata); 
		
		input = globaldata->ginput;
		list = globaldata->glist;
		order = list->getSharedOrderVector();
		
		while(order != NULL){
		
			if(globaldata->allLines == 1 || globaldata->lines.count(count) == 1 || globaldata->labels.count(order->getLabel()) == 1){
				//create collectors curve
				cCurve = new Collect(order, cDisplays);
				convert(globaldata->getFreq(), freq);
				cCurve->getSharedCurve(freq);
			
				delete cCurve;
			
				cout << order->getLabel() << '\t' << count << endl;
			}
			
			list = input->getListVector(); //get new list vector to process
			if (list != NULL) {
				order = list->getSharedOrderVector(); //gets new order vector with group info.
				count++;
			}else {
				break;
			}
		
		}
	
		for(int i=0;i<cDisplays.size();i++){	delete cDisplays[i];	}	
		return 0;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the CollectSharedCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the CollectSharedCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}


//**********************************************************************************************************************
