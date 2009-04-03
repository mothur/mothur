/*
 *  summarycommand.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "summarycommand.h"
#include "ace.h"
#include "sobs.h"
#include "nseqs.h"
#include "chao1.h"
#include "bootstrap.h"
#include "simpson.h"
#include "npshannon.h"
#include "shannon.h"
#include "jackknife.h"
#include "geom.h"
#include "logsd.h"
#include "qstat.h"
#include "bergerparker.h"
#include "bstick.h"

//**********************************************************************************************************************

SummaryCommand::SummaryCommand(){
	try {
		globaldata = GlobalData::getInstance();
		validCalculator = new ValidCalculators();
		int i;
		
		for (i=0; i<globaldata->Estimators.size(); i++) {
			if (validCalculator->isValidCalculator("summary", globaldata->Estimators[i]) == true) { 
				if(globaldata->Estimators[i] == "sobs"){
					sumCalculators.push_back(new Sobs());
				}else if(globaldata->Estimators[i] == "chao"){
					sumCalculators.push_back(new Chao1());
				}else if(globaldata->Estimators[i] == "geom"){
					sumCalculators.push_back(new Geom());
				}else if(globaldata->Estimators[i] == "logsd"){
					sumCalculators.push_back(new LogSD());
				}else if(globaldata->Estimators[i] == "qstat"){
					sumCalculators.push_back(new QStat());
				}else if(globaldata->Estimators[i] == "bergerparker"){
					sumCalculators.push_back(new BergerParker());
				}else if(globaldata->Estimators[i] == "bstick"){
					sumCalculators.push_back(new BStick());
				}else if(globaldata->Estimators[i] == "ace"){
					convert(globaldata->getAbund(), abund);
					if(abund < 5)
						abund = 10;
					sumCalculators.push_back(new Ace(abund));
				}else if(globaldata->Estimators[i] == "jack"){
					sumCalculators.push_back(new Jackknife());
				}else if(globaldata->Estimators[i] == "shannon"){
					sumCalculators.push_back(new Shannon());
				}else if(globaldata->Estimators[i] == "npshannon"){
					sumCalculators.push_back(new NPShannon());
				}else if(globaldata->Estimators[i] == "simpson"){
					sumCalculators.push_back(new Simpson());
				}else if(globaldata->Estimators[i] == "bootstrap"){
					sumCalculators.push_back(new Bootstrap());
				}else if (globaldata->Estimators[i] == "nseqs") { 
					sumCalculators.push_back(new NSeqs());
				}
			}
		}
		
		//reset calc for next command
		globaldata->setCalc("");

	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SummaryCommand class Function SummaryCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SummaryCommand class function SummaryCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}
//**********************************************************************************************************************

SummaryCommand::~SummaryCommand(){
	delete sabund;
	delete input;
	delete read;
}

//**********************************************************************************************************************

int SummaryCommand::execute(){
	try {
	
		//if the users entered no valid calculators don't execute command
		if (sumCalculators.size() == 0) { return 0; }

		outputFileName = ((getRootName(globaldata->inputFileName)) + "summary");
		openOutputFile(outputFileName, outputFileHandle);
		outputFileHandle << "label";
	
		read = new ReadPhilFile(globaldata->inputFileName);	
		read->read(&*globaldata); 

		for(int i=0;i<sumCalculators.size();i++){
			if(sumCalculators[i]->getCols() == 1){
				outputFileHandle << '\t' << sumCalculators[i]->getName();
			}
			else{
				outputFileHandle << '\t' << sumCalculators[i]->getName() << "\t" << sumCalculators[i]->getName() << "_lci\t" << sumCalculators[i]->getName() << "_hci";
			}
		}
		outputFileHandle << endl;
		
		sabund = globaldata->sabund;
		input = globaldata->ginput;
		int count = 1;
		while(sabund != NULL){
		
			if(globaldata->allLines == 1 || globaldata->lines.count(count) == 1 || globaldata->labels.count(sabund->getLabel()) == 1){			
	
				cout << sabund->getLabel() << '\t' << count << endl;
				
				outputFileHandle << sabund->getLabel();
				for(int i=0;i<sumCalculators.size();i++){
					vector<double> data = sumCalculators[i]->getValues(sabund);
					outputFileHandle << '\t';
					sumCalculators[i]->print(outputFileHandle);
				}
				
				outputFileHandle << endl;
			}
			sabund = input->getSAbundVector();
			count++;
		}
	
		return 0;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SummaryCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SummaryCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		
}

//**********************************************************************************************************************
