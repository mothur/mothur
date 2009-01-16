/*
 *  rarefactcommand.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "rarefactcommand.h"
#include "ace.h"
#include "sobs.h"
#include "chao1.h"
#include "bootstrap.h"
#include "simpson.h"
#include "npshannon.h"
#include "shannon.h"
#include "jackknife.h"

//**********************************************************************************************************************


RareFactCommand::RareFactCommand(){
	try {
		globaldata = GlobalData::getInstance();
		string fileNameRoot;
		fileNameRoot = getRootName(globaldata->inputFileName);
		int i;
		for (i=0; i<globaldata->rareEstimators.size(); i++) {
			if (globaldata->rareEstimators[i] == "sobs") { 
				rDisplays.push_back(new RareDisplay(new Sobs(), new ThreeColumnFile(fileNameRoot+"rarefaction")));
			}else if (globaldata->rareEstimators[i] == "chao") { 
				rDisplays.push_back(new RareDisplay(new Chao1(), new ThreeColumnFile(fileNameRoot+"r_chao")));
			}else if (globaldata->rareEstimators[i] == "ace") { 
				rDisplays.push_back(new RareDisplay(new Ace(), new ThreeColumnFile(fileNameRoot+"r_ace")));
			}else if (globaldata->rareEstimators[i] == "jack") { 
				rDisplays.push_back(new RareDisplay(new Jackknife(), new ThreeColumnFile(fileNameRoot+"r_jack")));
			}else if (globaldata->rareEstimators[i] == "shannon") { 
				rDisplays.push_back(new RareDisplay(new Shannon(), new ThreeColumnFile(fileNameRoot+"r_shannon")));
			}else if (globaldata->rareEstimators[i] == "npshannon") { 
				rDisplays.push_back(new RareDisplay(new NPShannon(), new ThreeColumnFile(fileNameRoot+"r_npshannon")));
			}else if (globaldata->rareEstimators[i] == "simpson") { 
				rDisplays.push_back(new RareDisplay(new Simpson(), new ThreeColumnFile(fileNameRoot+"r_simpson")));
			}else if (globaldata->rareEstimators[i] == "bootstrap") { 
				rDisplays.push_back(new RareDisplay(new Bootstrap(), new ThreeColumnFile(fileNameRoot+"r_bootstrap")));
			}
		}
	
	
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the RareFactCommand class Function RareFactCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the RareFactCommand class function RareFactCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
			
}

//**********************************************************************************************************************

RareFactCommand::~RareFactCommand(){
	delete order;
	delete input;
	delete rCurve;
	delete read;
}

//**********************************************************************************************************************

int RareFactCommand::execute(){
	try {
		int count = 1;
		read = new ReadPhilFile(globaldata->inputFileName);	
		read->read(&*globaldata); 

		order = globaldata->gorder;
		input = globaldata->ginput;
	
		while(order != NULL){
		
			if(globaldata->allLines == 1 || globaldata->lines.count(count) == 1 || globaldata->labels.count(order->getLabel()) == 1){
			
				rCurve = new Rarefact(order, rDisplays);
				convert(globaldata->getFreq(), freq);
				convert(globaldata->getIters(), nIters);
				rCurve->getCurve(freq, nIters);
			
				delete rCurve;
			
				cout << order->getLabel() << '\t' << count << endl;
			}
		
			order = (input->getOrderVector());
			count++;
		
		}
	
		for(int i=0;i<rDisplays.size();i++){	delete rDisplays[i];	}	
		return 0;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the RareFactCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the RareFactCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

//**********************************************************************************************************************
