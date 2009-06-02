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
	for (int i = 0; i < vennCalculators.size(); i++) {	delete vennCalculators[i];	}
}

//**********************************************************************************************************************

int VennCommand::execute(){
	try {
		int count = 1;
		SAbundVector* lastSAbund;
		vector<SharedRAbundVector*> lastLookup;	
		
		//if the users entered no valid calculators don't execute command
		if (vennCalculators.size() == 0) { return 0; }
		
		if (format == "sharedfile") {
			//you have groups
			read = new ReadOTUFile(globaldata->inputFileName);	
			read->read(&*globaldata); 
			
			input = globaldata->ginput;
			lookup = input->getSharedRAbundVectors();
			lastLookup = lookup;
		}else if (format == "list") {
			//you are using just a list file and have only one group
			read = new ReadOTUFile(globaldata->inputFileName);	
			read->read(&*globaldata); 
		
			sabund = globaldata->sabund;
			lastSAbund = globaldata->sabund;
			input = globaldata->ginput;
		}
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = globaldata->labels;
		
		if (format != "list") {	
			
			//as long as you are not at the end of the file or done wih the lines you want
			while((lookup[0] != NULL) && ((globaldata->allLines == 1) || (userLabels.size() != 0))) {
		
				if(globaldata->allLines == 1 || globaldata->lines.count(count) == 1 || globaldata->labels.count(lookup[0]->getLabel()) == 1){			
					cout << lookup[0]->getLabel() << '\t' << count << endl;
					processedLabels.insert(lookup[0]->getLabel());
					userLabels.erase(lookup[0]->getLabel());
					
					if (lookup.size() > 4) {
						cout << "Error: Too many groups chosen.  You may use up to 4 groups with the venn command.  I will use the first four groups in your groupfile." << endl;
						for (int i = lookup.size(); i > 4; i--) { lookup.pop_back(); } //no memmory leak because pop_back calls destructor
					}
					venn->getPic(lookup, vennCalculators);
				}
				
				if ((anyLabelsToProcess(lookup[0]->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLookup[0]->getLabel()) != 1)) {
					cout << lastLookup[0]->getLabel() << '\t' << count << endl;
					processedLabels.insert(lastLookup[0]->getLabel());
					userLabels.erase(lastLookup[0]->getLabel());

					if (lastLookup.size() > 4) {
						cout << "Error: Too many groups chosen.  You may use up to 4 groups with the venn command.  I will use the first four groups in your groupfile." << endl;
						for (int i = lastLookup.size(); i > 4; i--) { lastLookup.pop_back(); } //no memmory leak because pop_back calls destructor
					}				
					venn->getPic(lastLookup, vennCalculators);
				}
				
				//prevent memory leak
				if (count != 1) { for (int i = 0; i < lastLookup.size(); i++) {  delete lastLookup[i];  } }
				lastLookup = lookup;	
						
				//get next line to process
				lookup = input->getSharedRAbundVectors();
				count++;
			}
			
			//output error messages about any remaining user labels
			set<string>::iterator it;
			bool needToRun = false;
			for (it = userLabels.begin(); it != userLabels.end(); it++) {  
				cout << "Your file does not include the label "<< *it; 
				if (processedLabels.count(lastLookup[0]->getLabel()) != 1) {
					cout << ". I will use " << lastLookup[0]->getLabel() << "." << endl;
					needToRun = true;
				}else {
					cout << ". Please refer to " << lastLookup[0]->getLabel() << "." << endl;
				}
			}
		
			//run last line if you need to
			if (needToRun == true)  {
				cout << lastLookup[0]->getLabel() << '\t' << count << endl;
				if (lastLookup.size() > 4) {
					cout << "Error: Too many groups chosen.  You may use up to 4 groups with the venn command.  I will use the first four groups in your groupfile." << endl;
					for (int i = lastLookup.size(); i > 3; i--) { delete lastLookup[i]; lastLookup.pop_back(); }
				}				
				venn->getPic(lastLookup, vennCalculators);
			}
		
			for (int i = 0; i < lastLookup.size(); i++) {  delete lastLookup[i];  }

			//reset groups parameter
			globaldata->Groups.clear();  
			
		}else{
		
			while((sabund != NULL) && ((globaldata->allLines == 1) || (userLabels.size() != 0))) {
		
				if(globaldata->allLines == 1 || globaldata->lines.count(count) == 1 || globaldata->labels.count(sabund->getLabel()) == 1){			
	
					cout << sabund->getLabel() << '\t' << count << endl;
					venn->getPic(sabund, vennCalculators);
					
					processedLabels.insert(sabund->getLabel());
					userLabels.erase(sabund->getLabel());
				}
				
				if ((anyLabelsToProcess(sabund->getLabel(), userLabels, "") == true) && (processedLabels.count(lastSAbund->getLabel()) != 1)) {

					cout << lastSAbund->getLabel() << '\t' << count << endl;
					venn->getPic(lastSAbund, vennCalculators);
					
					processedLabels.insert(lastSAbund->getLabel());
					userLabels.erase(lastSAbund->getLabel());
				}		
				
				if (count != 1) { delete lastSAbund; }
				lastSAbund = sabund;			

				sabund = input->getSAbundVector();
				count++;
			}
			
			//output error messages about any remaining user labels
			set<string>::iterator it;
			bool needToRun = false;
			for (it = userLabels.begin(); it != userLabels.end(); it++) {  
				cout << "Your file does not include the label "<< *it; 
				if (processedLabels.count(lastSAbund->getLabel()) != 1) {
					cout << ". I will use " << lastSAbund->getLabel() << "." << endl;
					needToRun = true;
				}else {
					cout << ". Please refer to " << lastSAbund->getLabel() << "." << endl;
				}
			}
		
			//run last line if you need to
			if (needToRun == true)  {
				cout << lastSAbund->getLabel() << '\t' << count << endl;
				venn->getPic(lastSAbund, vennCalculators);
			}
			delete lastSAbund;
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