/*
 *  summarysharedcommand.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "summarysharedcommand.h"
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

SummarySharedCommand::SummarySharedCommand(){
	try {
		globaldata = GlobalData::getInstance();
		
		int i;
		for (i=0; i<globaldata->sharedSummaryEstimators.size(); i++) {
			if (globaldata->sharedSummaryEstimators[i] == "sharedChao") { 
				sumCalculators.push_back(new SharedChao1());
			}else if (globaldata->sharedSummaryEstimators[i] == "sharedAce") { 
				sumCalculators.push_back(new SharedAce());
			}else if (globaldata->sharedSummaryEstimators[i] == "sharedJabund") { 	
				sumCalculators.push_back(new SharedJAbund());
			}else if (globaldata->sharedSummaryEstimators[i] == "sharedSorensonAbund") { 
				sumCalculators.push_back(new SharedSorAbund());
			}else if (globaldata->sharedSummaryEstimators[i] == "sharedJclass") { 
				sumCalculators.push_back(new SharedJclass());
			}else if (globaldata->sharedSummaryEstimators[i] == "sharedSorClass") { 
				sumCalculators.push_back(new SharedSorClass());
			}else if (globaldata->sharedSummaryEstimators[i] == "sharedJest") { 
				sumCalculators.push_back(new SharedJest());
			}else if (globaldata->sharedSummaryEstimators[i] == "sharedSorEst") { 
				sumCalculators.push_back(new SharedSorEst());
			}else if (globaldata->sharedSummaryEstimators[i] == "SharedThetaYC") { 
				sumCalculators.push_back(new SharedThetaYC());
			}else if (globaldata->sharedSummaryEstimators[i] == "SharedThetaN") { 
				sumCalculators.push_back(new SharedThetaN());
			}
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SummarySharedCommand class Function SummarySharedCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SummarySharedCommand class function SummarySharedCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}
//**********************************************************************************************************************

SummarySharedCommand::~SummarySharedCommand(){
	delete input;
	delete read;
}

//**********************************************************************************************************************

int SummarySharedCommand::execute(){
	try {
		outputFileName = ((getRootName(globaldata->inputFileName)) + "sharedSummary");
		openOutputFile(outputFileName, outputFileHandle);
	
		read = new ReadPhilFile(globaldata->inputFileName);	
		read->read(&*globaldata); 
		
		outputFileHandle << '\t' << '\t' << '\t' << '\t'; //pads file for labels and groupnames
		for(int i=0;i<sumCalculators.size();i++){
			outputFileHandle << '\t' << sumCalculators[i]->getName();
		}
		outputFileHandle << endl;
		
		list = globaldata->glist;
		input = globaldata->ginput;
		order = list->getSharedOrderVector();
		getGroupComb();
		
		int count = 1;
		while(order != NULL){
		
			if(globaldata->allLines == 1 || globaldata->lines.count(count) == 1 || globaldata->labels.count(order->getLabel()) == 1){			
	
				cout << order->getLabel() << '\t' << count << endl;
				getSharedVectors();  //fills group vectors from order vector.
				
				//randomize group order
				if (globaldata->getJumble() == "1") { random_shuffle(lookup.begin(), lookup.end()); }

				int n = 1; 
				for (int k = 0; k < (lookup.size() - 1); k++) { // pass cdd each set of groups to commpare
					for (int l = n; l < lookup.size(); l++) {
						outputFileHandle << order->getLabel() << '\t' << groupComb[n-1] << '\t'; //print out label and group
						for(int i=0;i<sumCalculators.size();i++){
							sumCalculators[i]->getValues(lookup[k], lookup[l]); //saves the calculator outputs
							outputFileHandle << '\t';
							sumCalculators[i]->print(outputFileHandle);
						}
						outputFileHandle << endl;
					}
					n++;
				}
			}
		
			list = input->getListVector(); //get new list vector to process
			if (list != NULL) {
				order = list->getSharedOrderVector(); //gets new order vector with group info.
				count++;
			}else {
				break;
			}
		}
	
		return 0;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SummarySharedCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SummarySharedCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		
}

//**********************************************************************************************************************

void SummarySharedCommand::getSharedVectors(){
try {
		lookup.clear();
		//create and initialize vector of sharedvectors, one for each group
		for (int i = 0; i < globaldata->gGroupmap->getNumGroups(); i++) { 
			SharedRAbundVector* temp = new SharedRAbundVector(order->getNumBins());
			temp->setLabel(order->getLabel());
			temp->setGroup(globaldata->gGroupmap->namesOfGroups[i]);
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
		cout << "Standard Error: " << e.what() << " has occurred in the SummarySharedCommand class Function getSharedVectors. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SummarySharedCommand class function getSharedVectors. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}

}

/**************************************************************************************/
void SummarySharedCommand::getGroupComb() {
	try {
		string group;
		
		int n = 1;
		for (int i = 0; i < (globaldata->gGroupmap->getNumGroups() - 1); i++) {
			for (int l = n; l < globaldata->gGroupmap->getNumGroups(); l++) {
				group = globaldata->gGroupmap->namesOfGroups[i] + globaldata->gGroupmap->namesOfGroups[l];
				groupComb.push_back(group);	
			}
			n++;
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SummarySharedCommand class Function getGroupComb. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SummarySharedCommand class function getGroupComb. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}

}

