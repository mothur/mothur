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
#include "whittaker.h"
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
		mult = false;
		
		int i;
		for (i=0; i<globaldata->Estimators.size(); i++) {
			if (validCalculator->isValidCalculator("sharedsummary", globaldata->Estimators[i]) == true) { 
				if (globaldata->Estimators[i] == "sharedsobs") { 
					sumCalculators.push_back(new SharedSobsCS());
				}else if (globaldata->Estimators[i] == "sharedchao") { 
					sumCalculators.push_back(new SharedChao1());
				}else if (globaldata->Estimators[i] == "sharedace") { 
					sumCalculators.push_back(new SharedAce());
				}else if (globaldata->Estimators[i] == "jabund") { 	
					sumCalculators.push_back(new JAbund());
				}else if (globaldata->Estimators[i] == "sorabund") { 
					sumCalculators.push_back(new SorAbund());
				}else if (globaldata->Estimators[i] == "jclass") { 
					sumCalculators.push_back(new Jclass());
				}else if (globaldata->Estimators[i] == "sorclass") { 
					sumCalculators.push_back(new SorClass());
				}else if (globaldata->Estimators[i] == "jest") { 
					sumCalculators.push_back(new Jest());
				}else if (globaldata->Estimators[i] == "sorest") { 
					sumCalculators.push_back(new SorEst());
				}else if (globaldata->Estimators[i] == "thetayc") { 
					sumCalculators.push_back(new ThetaYC());
				}else if (globaldata->Estimators[i] == "thetan") { 
					sumCalculators.push_back(new ThetaN());
				}else if (globaldata->Estimators[i] == "kstest") { 
					sumCalculators.push_back(new KSTest());
				}else if (globaldata->Estimators[i] == "sharednseqs") { 
					sumCalculators.push_back(new SharedNSeqs());
				}else if (globaldata->Estimators[i] == "ochiai") { 
					sumCalculators.push_back(new Ochiai());
				}else if (globaldata->Estimators[i] == "anderberg") { 
					sumCalculators.push_back(new Anderberg());
				}else if (globaldata->Estimators[i] == "kulczynski") { 
					sumCalculators.push_back(new Kulczynski());
				}else if (globaldata->Estimators[i] == "kulczynskicody") { 
					sumCalculators.push_back(new KulczynskiCody());
				}else if (globaldata->Estimators[i] == "lennon") { 
					sumCalculators.push_back(new Lennon());
				}else if (globaldata->Estimators[i] == "morisitahorn") { 
					sumCalculators.push_back(new MorHorn());
				}else if (globaldata->Estimators[i] == "braycurtis") { 
					sumCalculators.push_back(new BrayCurtis());
				}else if (globaldata->Estimators[i] == "whittaker") { 
					sumCalculators.push_back(new Whittaker());
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
		//check if any calcs can do multiples
		else{
			for (int i = 0; i < sumCalculators.size(); i++) {
				if (sumCalculators[i]->getMultiple() == true) { mult = true; }
			}
		}
		
		if (format == "sharedfile") {
			read = new ReadOTUFile(globaldata->inputFileName);	
			read->read(&*globaldata); 
			
			input = globaldata->ginput;
			order = input->getSharedOrderVector();
		}else {
			//you are using a list and a groupfile
			read = new ReadOTUFile(globaldata->inputFileName);	
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
		
		//create file and put column headers for multiple groups file
		if (mult = true) {
			outAllFileName = ((getRootName(globaldata->inputFileName)) + "sharedmultiple.summary");
			openOutputFile(outAllFileName, outAll);
			
			outAll << "label" <<'\t' << "comparison" << '\t'; 
			for(int i=0;i<sumCalculators.size();i++){
				if (sumCalculators[i]->getMultiple() == true) { 
					outAll << '\t' << sumCalculators[i]->getName();
				}
			}
			outAll << endl;
		}
		
		while(order != NULL){
		
			if(globaldata->allLines == 1 || globaldata->lines.count(count) == 1 || globaldata->labels.count(order->getLabel()) == 1){			
	
				cout << order->getLabel() << '\t' << count << endl;
				util->getSharedVectors(globaldata->Groups, lookup, order);  //fills group vectors from order vector.  //fills group vectors from order vector.
				
				//loop through calculators and add to file all for all calcs that can do mutiple groups
				if (mult = true) {
					//output label
					outAll << order->getLabel() << '\t';
					
					//output groups names
					string outNames = "";
					for (int j = 0; j < lookup.size(); j++) {
						outNames += lookup[j]->getGroup() +  "-";
					}
					outNames = outNames.substr(0, outNames.length()-1); //rip off extra '-';
					outAll << outNames << '\t';
					
					for(int i=0;i<sumCalculators.size();i++){
						if (sumCalculators[i]->getMultiple() == true) { 
							sumCalculators[i]->getValues(lookup);
							outAll << '\t';
							sumCalculators[i]->print(outAll);
						}
					}
					outAll << endl;
				}

				int n = 1; 
				vector<SharedRAbundVector*> subset;
				for (int k = 0; k < (lookup.size() - 1); k++) { // pass cdd each set of groups to commpare
					for (int l = n; l < lookup.size(); l++) {
						
						outputFileHandle << order->getLabel() << '\t';
						
						subset.clear(); //clear out old pair of sharedrabunds
						//add new pair of sharedrabunds
						subset.push_back(lookup[k]); subset.push_back(lookup[l]); 
						
						//sort groups to be alphanumeric
						if (lookup[k]->getGroup() > lookup[l]->getGroup()) {
							outputFileHandle << (lookup[l]->getGroup() +'\t' + lookup[k]->getGroup()) << '\t'; //print out groups
						}else{
							outputFileHandle << (lookup[k]->getGroup() +'\t' + lookup[l]->getGroup()) << '\t'; //print out groups
						}
						
						for(int i=0;i<sumCalculators.size();i++) {
							sumCalculators[i]->getValues(subset); //saves the calculator outputs
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
