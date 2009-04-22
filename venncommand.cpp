/*
 *  venncommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 3/30/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "venncommand.h"
#include "ace.h"
#include "sobs.h"
#include "chao1.h"
//#include "jackknife.h"
#include "sharedsobscollectsummary.h"
#include "sharedchao1.h"
#include "sharedace.h"


//**********************************************************************************************************************

VennCommand::VennCommand(){
	try {
		globaldata = GlobalData::getInstance();
		format = globaldata->getFormat();
		validCalculator = new ValidCalculators();
		
		int i;
		
		if (format == "list") {
			for (i=0; i<globaldata->Estimators.size(); i++) {
				if (validCalculator->isValidCalculator("vennsingle", globaldata->Estimators[i]) == true) { 
					if (globaldata->Estimators[i] == "sobs") { 
						vennCalculators.push_back(new Sobs());
					}else if (globaldata->Estimators[i] == "chao") { 
						vennCalculators.push_back(new Chao1());
					}else if (globaldata->Estimators[i] == "ace") {
						convert(globaldata->getAbund(), abund);
						if(abund < 5)
							abund = 10;
						vennCalculators.push_back(new Ace(abund));
					//}else if (globaldata->Estimators[i] == "jack") { 	
						//vennCalculators.push_back(new Jackknife());
					}
				}
			}
		}else {
			for (i=0; i<globaldata->Estimators.size(); i++) {
				if (validCalculator->isValidCalculator("vennshared", globaldata->Estimators[i]) == true) { 
					if (globaldata->Estimators[i] == "sharedsobs") { 
						vennCalculators.push_back(new SharedSobsCS());
					}else if (globaldata->Estimators[i] == "sharedchao") { 
						vennCalculators.push_back(new SharedChao1());
					}else if (globaldata->Estimators[i] == "sharedace") { 
						vennCalculators.push_back(new SharedAce());
					}
				}
			}
		}
		
		venn = new Venn();
		
		//reset calc for next command
		globaldata->setCalc("");

		
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the VennCommand class Function VennCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the VennCommand class function VennCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}
//**********************************************************************************************************************

VennCommand::~VennCommand(){
	delete input;
	delete read;
	delete venn;
}

//**********************************************************************************************************************

int VennCommand::execute(){
	try {
		int count = 1;	
		
		//if the users entered no valid calculators don't execute command
		if (vennCalculators.size() == 0) { return 0; }
		
		if (format == "sharedfile") {
			//you have groups
			read = new ReadOTUFile(globaldata->inputFileName);	
			read->read(&*globaldata); 
			
			input = globaldata->ginput;
			lookup = input->getSharedRAbundVectors();
		}else if (format == "shared") {
			//you are using a list and a groupfile
			read = new ReadOTUFile(globaldata->inputFileName);	
			read->read(&*globaldata); 
		
			input = globaldata->ginput;
			SharedList = globaldata->gSharedList;
			lookup = SharedList->getSharedRAbundVector();
		}else if (format == "list") {
			//you are using just a list file and have only one group
			read = new ReadOTUFile(globaldata->inputFileName);	
			read->read(&*globaldata); 
		
			sabund = globaldata->sabund;
			input = globaldata->ginput;
		}

		
		if (format != "list") {	
			
			while(lookup[0] != NULL){
		
				if(globaldata->allLines == 1 || globaldata->lines.count(count) == 1 || globaldata->labels.count(lookup[0]->getLabel()) == 1){			
	
					cout << lookup[0]->getLabel() << '\t' << count << endl;
					
					if (lookup.size() > 4) {
						cout << "Error: Too many groups chosen.  You may use up to 4 groups with the venn command.  I will use the first four groups in your groupfile." << endl;
						for (int i = lookup.size(); i > 3; i--) { delete lookup[i]; lookup.pop_back(); }
					}
					
					//util->getSharedVectors(globaldata->Groups, lookup, order);  //fills group vectors from order vector.
					venn->getPic(lookup, vennCalculators);
				}
						
				//get next line to process
				lookup = input->getSharedRAbundVectors();
				count++;
			}
			
			//reset groups parameter
			globaldata->Groups.clear();  
			
		}else{
		
			while(sabund != NULL){
		
				if(globaldata->allLines == 1 || globaldata->lines.count(count) == 1 || globaldata->labels.count(sabund->getLabel()) == 1){			
	
					cout << sabund->getLabel() << '\t' << count << endl;
					venn->getPic(sabund, vennCalculators);
				}
				
				sabund = input->getSAbundVector();
				count++;
			}
		}
		
		globaldata->setGroups("");
		return 0;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the VennCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the VennCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		
}

//**********************************************************************************************************************