/*
 *  summarysharedcommand.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "summarysharedcommand.h"
#include "sharedsobscollectsummary.h"
#include "sharedchao1.h"
#include "sharedace.h"
#include "sharednseqs.h"
#include "sharedjabund.h"
#include "sharedsorabund.h"
#include "sharedjclass.h"
#include "sharedsorclass.h"
#include "sharedjest.h"
#include "sharedsorest.h"
#include "sharedthetayc.h"
#include "sharedthetan.h"
#include "sharedkstest.h"
#include "sharedbdiversity.h"
#include "sharedochiai.h"
#include "sharedanderbergs.h"
#include "sharedkulczynski.h"
#include "sharedkulczynskicody.h"
#include "sharedlennon.h"
#include "sharedmorisitahorn.h"
#include "sharedbraycurtis.h"


//**********************************************************************************************************************

SummarySharedCommand::SummarySharedCommand(){
	try {
		globaldata = GlobalData::getInstance();
		outputFileName = ((getRootName(globaldata->inputFileName)) + "shared.summary");
		openOutputFile(outputFileName, outputFileHandle);
		format = globaldata->getFormat();
		validCalculator = new ValidCalculators();
		util = new SharedUtil();
		
		int i;
		for (i=0; i<globaldata->Estimators.size(); i++) {
			if (validCalculator->isValidCalculator("sharedsummary", globaldata->Estimators[i]) == true) { 
				if (globaldata->Estimators[i] == "sharedsobs") { 
					sumCalculators.push_back(new SharedSobsCS());
				}else if (globaldata->Estimators[i] == "sharedchao") { 
					sumCalculators.push_back(new SharedChao1());
				}else if (globaldata->Estimators[i] == "sharedace") { 
					sumCalculators.push_back(new SharedAce());
				}else if (globaldata->Estimators[i] == "sharedjabund") { 	
					sumCalculators.push_back(new SharedJAbund());
				}else if (globaldata->Estimators[i] == "sharedsorensonabund") { 
					sumCalculators.push_back(new SharedSorAbund());
				}else if (globaldata->Estimators[i] == "sharedjclass") { 
					sumCalculators.push_back(new SharedJclass());
				}else if (globaldata->Estimators[i] == "sharedsorclass") { 
					sumCalculators.push_back(new SharedSorClass());
				}else if (globaldata->Estimators[i] == "sharedjest") { 
					sumCalculators.push_back(new SharedJest());
				}else if (globaldata->Estimators[i] == "sharedsorest") { 
					sumCalculators.push_back(new SharedSorEst());
				}else if (globaldata->Estimators[i] == "sharedthetayc") { 
					sumCalculators.push_back(new SharedThetaYC());
				}else if (globaldata->Estimators[i] == "sharedthetan") { 
					sumCalculators.push_back(new SharedThetaN());
				}else if (globaldata->Estimators[i] == "sharedkstest") { 
					sumCalculators.push_back(new SharedKSTest());
				}else if (globaldata->Estimators[i] == "sharednseqs") { 
					sumCalculators.push_back(new SharedNSeqs());
				}else if (globaldata->Estimators[i] == "sharedochiai") { 
					sumCalculators.push_back(new SharedOchiai());
				}else if (globaldata->Estimators[i] == "sharedanderberg") { 
					sumCalculators.push_back(new SharedAnderberg());
				}else if (globaldata->Estimators[i] == "sharedkulczynski") { 
					sumCalculators.push_back(new SharedKulczynski());
				}else if (globaldata->Estimators[i] == "sharedkulczynskicody") { 
					sumCalculators.push_back(new SharedKulczynskiCody());
				}else if (globaldata->Estimators[i] == "sharedlennon") { 
					sumCalculators.push_back(new SharedLennon());
				}else if (globaldata->Estimators[i] == "sharedmorisitahorn") { 
					sumCalculators.push_back(new SharedMorHorn());
				}else if (globaldata->Estimators[i] == "sharedbraycurtis") { 
					sumCalculators.push_back(new SharedBrayCurtis());
				}
				else if (globaldata->Estimators[i] == "sharedbdiversity") { 
					sumCalculators.push_back(new SharedBDiversity());
				}
				
			}
		}
		//reset calc for next command
		globaldata->setCalc("");

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
	delete util;
}

//**********************************************************************************************************************

int SummarySharedCommand::execute(){
	try {
		int count = 1;	
		
		//if the users entered no valid calculators don't execute command
		if (sumCalculators.size() == 0) { return 0; }

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
		util->setGroups(globaldata->Groups, globaldata->gGroupmap->namesOfGroups, "summary");
		
		//output estimator names as column headers
		outputFileHandle << "label" <<'\t' << "comparison" << '\t'; 
		for(int i=0;i<sumCalculators.size();i++){
			outputFileHandle << '\t' << sumCalculators[i]->getName();
		}
		outputFileHandle << endl;
		
		while(order != NULL){
		
			if(globaldata->allLines == 1 || globaldata->lines.count(count) == 1 || globaldata->labels.count(order->getLabel()) == 1){			
	
				cout << order->getLabel() << '\t' << count << endl;
				util->getSharedVectors(globaldata->Groups, lookup, order);  //fills group vectors from order vector.  //fills group vectors from order vector.
				
				//randomize group order
				if (globaldata->getJumble() == "1") { random_shuffle(lookup.begin(), lookup.end()); }

				int n = 1; 
				for (int k = 0; k < (lookup.size() - 1); k++) { // pass cdd each set of groups to commpare
					for (int l = n; l < lookup.size(); l++) {
						
						outputFileHandle << order->getLabel() << '\t';
						
						//sort groups to be alphanumeric
						if (lookup[k]->getGroup() > lookup[l]->getGroup()) {
							outputFileHandle << (lookup[l]->getGroup() +'\t' + lookup[k]->getGroup()) << '\t'; //print out groups
						}else{
							outputFileHandle << (lookup[k]->getGroup() +'\t' + lookup[l]->getGroup()) << '\t'; //print out groups
						}
						
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

/***********************************************************/
