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
#include "sharedjackknife.h"
#include "whittaker.h"


//**********************************************************************************************************************

SummarySharedCommand::SummarySharedCommand(){
	try {
		globaldata = GlobalData::getInstance();
		outputFileName = ((getRootName(globaldata->inputFileName)) + "shared.summary");
		openOutputFile(outputFileName, outputFileHandle);
		format = globaldata->getFormat();
		validCalculator = new ValidCalculators();
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
		
		//read first line
		read = new ReadOTUFile(globaldata->inputFileName);	
		read->read(&*globaldata); 
			
		input = globaldata->ginput;
		lookup = input->getSharedRAbundVectors();
		vector<SharedRAbundVector*> lastLookup = lookup;
		
		//output estimator names as column headers
		outputFileHandle << "label" <<'\t' << "comparison" << '\t'; 
		for(int i=0;i<sumCalculators.size();i++){
			outputFileHandle << '\t' << sumCalculators[i]->getName();
		}
		outputFileHandle << endl;
		
		//create file and put column headers for multiple groups file
		if (mult == true) {
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
		
		if (lookup.size() < 2) { 
			cout << "I cannot run the command without at least 2 valid groups."; 
			for (int i = 0; i < lookup.size(); i++) { delete lookup[i]; }
			
			//close files and clean up
			outputFileHandle.close();  remove(outputFileName.c_str());
			if (mult == true) {  outAll.close();  remove(outAllFileName.c_str());  }
			return 0;
		}
					
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = globaldata->labels;
		set<int> userLines = globaldata->lines;
		
		//as long as you are not at the end of the file or done wih the lines you want
		while((lookup[0] != NULL) && ((globaldata->allLines == 1) || (userLabels.size() != 0) || (userLines.size() != 0))) {
		
			if(globaldata->allLines == 1 || globaldata->lines.count(count) == 1 || globaldata->labels.count(lookup[0]->getLabel()) == 1){			
				cout << lookup[0]->getLabel() << '\t' << count << endl;
				process(lookup);
				
				processedLabels.insert(lookup[0]->getLabel());
				userLabels.erase(lookup[0]->getLabel());
				userLines.erase(count);
			}
			
			if ((anyLabelsToProcess(lookup[0]->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLookup[0]->getLabel()) != 1)) {
					cout << lastLookup[0]->getLabel() << '\t' << count << endl;
					process(lastLookup);
					
					processedLabels.insert(lastLookup[0]->getLabel());
					userLabels.erase(lastLookup[0]->getLabel());
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
			process(lastLookup);
		}
		
		for (int i = 0; i < lastLookup.size(); i++) {  delete lastLookup[i];  }

		//reset groups parameter
		globaldata->Groups.clear();  globaldata->setGroups("");
		
		//close files
		outputFileHandle.close();
		if (mult == true) {  outAll.close();  }

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
void SummarySharedCommand::process(vector<SharedRAbundVector*> thisLookup) {
	try {
				//loop through calculators and add to file all for all calcs that can do mutiple groups
				if (mult == true) {
					//output label
					outAll << thisLookup[0]->getLabel() << '\t';
					
					//output groups names
					string outNames = "";
					for (int j = 0; j < thisLookup.size(); j++) {
						outNames += thisLookup[j]->getGroup() +  "-";
					}
					outNames = outNames.substr(0, outNames.length()-1); //rip off extra '-';
					outAll << outNames << '\t';
					
					for(int i=0;i<sumCalculators.size();i++){
						if (sumCalculators[i]->getMultiple() == true) { 
							sumCalculators[i]->getValues(thisLookup);
							outAll << '\t';
							sumCalculators[i]->print(outAll);
						}
					}
					outAll << endl;
				}
	
				int n = 1; 
				vector<SharedRAbundVector*> subset;
				for (int k = 0; k < (thisLookup.size() - 1); k++) { // pass cdd each set of groups to commpare
					for (int l = n; l < thisLookup.size(); l++) {
						
						outputFileHandle << thisLookup[0]->getLabel() << '\t';
						
						subset.clear(); //clear out old pair of sharedrabunds
						//add new pair of sharedrabunds
						subset.push_back(thisLookup[k]); subset.push_back(thisLookup[l]); 
						
						//sort groups to be alphanumeric
						if (thisLookup[k]->getGroup() > thisLookup[l]->getGroup()) {
							outputFileHandle << (thisLookup[l]->getGroup() +'\t' + thisLookup[k]->getGroup()) << '\t'; //print out groups
						}else{
							outputFileHandle << (thisLookup[k]->getGroup() +'\t' + thisLookup[l]->getGroup()) << '\t'; //print out groups
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
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SummarySharedCommand class Function process. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SummarySharedCommand class function process. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		
}

/***********************************************************/