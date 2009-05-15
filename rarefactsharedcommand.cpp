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
#include "sharednseqs.h"

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
				}else if (globaldata->Estimators[i] == "sharednseqs") { 
					rDisplays.push_back(new RareDisplay(new SharedNSeqs(), new SharedThreeColumnFile(fileNameRoot+"shared.r_nseqs", "")));
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

		read = new ReadOTUFile(globaldata->inputFileName);	
		read->read(&*globaldata); 
			
		input = globaldata->ginput;
		lookup = input->getSharedRAbundVectors();
		
		if (lookup.size() < 2) { 
			cout << "I cannot run the command without at least 2 valid groups."; 
			for (int i = 0; i < lookup.size(); i++) { delete lookup[i]; }
			return 0;
		}
					
		
		while(lookup[0] != NULL){
		
			if(globaldata->allLines == 1 || globaldata->lines.count(count) == 1 || globaldata->labels.count(lookup[0]->getLabel()) == 1){
				//create collectors curve
				rCurve = new Rarefact(lookup, rDisplays);
				convert(globaldata->getFreq(), freq);
				convert(globaldata->getIters(), nIters);
				rCurve->getSharedCurve(freq, nIters);
			
				delete rCurve;
			
				cout << lookup[0]->getLabel() << '\t' << count << endl;
			}
			
			//prevent memory leak
			for (int i = 0; i < lookup.size(); i++) {	delete lookup[i];	}
				
			//get next line to process
			lookup = input->getSharedRAbundVectors();
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
