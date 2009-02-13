/*
 *  collectcommand.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "collectcommand.h"
#include "ace.h"
#include "sobs.h"
#include "chao1.h"
#include "bootstrap.h"
#include "simpson.h"
#include "npshannon.h"
#include "shannon.h"
#include "jackknife.h"


//**********************************************************************************************************************


CollectCommand::CollectCommand(){
	try {
		globaldata = GlobalData::getInstance();
		string fileNameRoot;
		fileNameRoot = getRootName(globaldata->inputFileName);
		int i;
		validCalculator = new ValidCalculators();
		
		for (i=0; i<globaldata->Estimators.size(); i++) {
			if (validCalculator->isValidCalculator("single", globaldata->Estimators[i]) == true) { 
				if (globaldata->Estimators[i] == "sobs") { 
					cDisplays.push_back(new CollectDisplay(new Sobs(), new OneColumnFile(fileNameRoot+"sobs")));
				}else if (globaldata->Estimators[i] == "chao") { 
					cDisplays.push_back(new CollectDisplay(new Chao1(), new ThreeColumnFile(fileNameRoot+"chao")));
				}else if (globaldata->Estimators[i] == "ace") { 
					cDisplays.push_back(new CollectDisplay(new Ace(), new ThreeColumnFile(fileNameRoot+"ace")));
				}else if (globaldata->Estimators[i] == "jack") { 
					cDisplays.push_back(new CollectDisplay(new Jackknife(), new ThreeColumnFile(fileNameRoot+"jack")));
				}else if (globaldata->Estimators[i] == "shannon") { 
					cDisplays.push_back(new CollectDisplay(new Shannon(), new ThreeColumnFile(fileNameRoot+"shannon")));
				}else if (globaldata->Estimators[i] == "npshannon") { 
					cDisplays.push_back(new CollectDisplay(new NPShannon(), new OneColumnFile(fileNameRoot+"np_shannon")));
				}else if (globaldata->Estimators[i] == "simpson") { 
					cDisplays.push_back(new CollectDisplay(new Simpson(), new ThreeColumnFile(fileNameRoot+"simpson")));
				}else if (globaldata->Estimators[i] == "bootstrap") { 
					cDisplays.push_back(new CollectDisplay(new Bootstrap(), new OneColumnFile(fileNameRoot+"bootstrap")));
				}
			}
		}
		
		//reset calc for next command
		globaldata->setCalc("");
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the CollectCommand class Function CollectCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the CollectCommand class function CollectCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
			
}

//**********************************************************************************************************************

CollectCommand::~CollectCommand(){
	delete order;
	delete input;
	delete cCurve;
	delete read;
}

//**********************************************************************************************************************

int CollectCommand::execute(){
	try {
		int count = 1;
		
		//if the users entered no valid calculators don't execute command
		if (cDisplays.size() == 0) { return 0; }

		read = new ReadPhilFile(globaldata->inputFileName);	
		read->read(&*globaldata); 
		
		order = globaldata->gorder;
		input = globaldata->ginput;
			
		while(order != NULL){
		
			if(globaldata->allLines == 1 || globaldata->lines.count(count) == 1 || globaldata->labels.count(order->getLabel()) == 1){
			
				cCurve = new Collect(order, cDisplays);
				convert(globaldata->getFreq(), freq);
				cCurve->getCurve(freq);
			
				delete cCurve;
			
				cout << order->getLabel() << '\t' << count << endl;
			}
		
			order = (input->getOrderVector());
			count++;
		
		}
	
		for(int i=0;i<cDisplays.size();i++){	delete cDisplays[i];	}	
		return 0;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the CollectCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the CollectCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

//**********************************************************************************************************************
