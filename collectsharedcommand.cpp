/*
 *  collectsharedcommand.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "collectsharedcommand.h"
#include "sharedsobscollectsummary.h"
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
		format = globaldata->getFormat();
		validCalculator = new ValidCalculators();
		
		int i;
		for (i=0; i<globaldata->Estimators.size(); i++) {
			if (validCalculator->isValidCalculator("shared", globaldata->Estimators[i]) == true) { 
				if (globaldata->Estimators[i] == "sharedchao") { 
					cDisplays.push_back(new CollectDisplay(new SharedChao1(), new SharedOneColumnFile(fileNameRoot+"shared.chao")));
				}else if (globaldata->Estimators[i] == "sharedsobs") { 
					cDisplays.push_back(new CollectDisplay(new SharedSobsCS(), new SharedOneColumnFile(fileNameRoot+"shared.sobs")));
				}else if (globaldata->Estimators[i] == "sharedace") { 
					cDisplays.push_back(new CollectDisplay(new SharedAce(), new SharedOneColumnFile(fileNameRoot+"shared.ace")));
				}else if (globaldata->Estimators[i] == "sharedjabund") { 	
					cDisplays.push_back(new CollectDisplay(new SharedJAbund(), new SharedOneColumnFile(fileNameRoot+"shared.jabund")));
				}else if (globaldata->Estimators[i] == "sharedsorensonabund") { 
					cDisplays.push_back(new CollectDisplay(new SharedSorAbund(), new SharedOneColumnFile(fileNameRoot+"shared.sorabund")));
				}else if (globaldata->Estimators[i] == "sharedjclass") { 
					cDisplays.push_back(new CollectDisplay(new SharedJclass(), new SharedOneColumnFile(fileNameRoot+"shared.jclass")));
				}else if (globaldata->Estimators[i] == "sharedsorclass") { 
					cDisplays.push_back(new CollectDisplay(new SharedSorClass(), new SharedOneColumnFile(fileNameRoot+"shared.sorclass")));
				}else if (globaldata->Estimators[i] == "sharedjest") { 
					cDisplays.push_back(new CollectDisplay(new SharedJest(), new SharedOneColumnFile(fileNameRoot+"shared.jest")));
				}else if (globaldata->Estimators[i] == "sharedsorest") { 
					cDisplays.push_back(new CollectDisplay(new SharedSorEst(), new SharedOneColumnFile(fileNameRoot+"shared.sorest")));
				}else if (globaldata->Estimators[i] == "sharedthetayc") { 
					cDisplays.push_back(new CollectDisplay(new SharedThetaYC(), new SharedOneColumnFile(fileNameRoot+"shared.thetayc")));
				}else if (globaldata->Estimators[i] == "sharedthetan") { 
					cDisplays.push_back(new CollectDisplay(new SharedThetaN(), new SharedOneColumnFile(fileNameRoot+"shared.thetan")));
				}
			}
		}
		
		//reset calc for next command
		globaldata->setCalc("");

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
		
		//if the users entered no valid calculators don't execute command
		if (cDisplays.size() == 0) { return 0; }
		
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
		
		while(order != NULL){
		
			if(globaldata->allLines == 1 || globaldata->lines.count(count) == 1 || globaldata->labels.count(order->getLabel()) == 1){
				//create collectors curve
				cCurve = new Collect(order, cDisplays);
				convert(globaldata->getFreq(), freq);
				cCurve->getSharedCurve(freq);
			
				delete cCurve;
			
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
