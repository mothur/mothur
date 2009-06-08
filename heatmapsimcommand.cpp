/*
 *  heatmapsimcommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 6/8/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "heatmapsimcommand.h"
#include "sharedjabund.h"
#include "sharedsorabund.h"
#include "sharedjclass.h"
#include "sharedsorclass.h"
#include "sharedjest.h"
#include "sharedsorest.h"
#include "sharedthetayc.h"
#include "sharedthetan.h"
#include "sharedmorisitahorn.h"
#include "sharedbraycurtis.h"


//**********************************************************************************************************************

HeatMapSimCommand::HeatMapSimCommand(){
	try {
		globaldata = GlobalData::getInstance();
		validCalculator = new ValidCalculators();
		heatmap = new HeatMapSim();
			
		int i;
		for (i=0; i<globaldata->Estimators.size(); i++) {
			if (validCalculator->isValidCalculator("heat", globaldata->Estimators[i]) == true) { 
				if (globaldata->Estimators[i] == "jabund") { 	
					heatCalculators.push_back(new JAbund());
				}else if (globaldata->Estimators[i] == "sorabund") { 
					heatCalculators.push_back(new SorAbund());
				}else if (globaldata->Estimators[i] == "jclass") { 
					heatCalculators.push_back(new Jclass());
				}else if (globaldata->Estimators[i] == "sorclass") { 
					heatCalculators.push_back(new SorClass());
				}else if (globaldata->Estimators[i] == "jest") { 
					heatCalculators.push_back(new Jest());
				}else if (globaldata->Estimators[i] == "sorest") { 
					heatCalculators.push_back(new SorEst());
				}else if (globaldata->Estimators[i] == "thetayc") { 
					heatCalculators.push_back(new ThetaYC());
				}else if (globaldata->Estimators[i] == "thetan") { 
					heatCalculators.push_back(new ThetaN());
				}else if (globaldata->Estimators[i] == "morisitahorn") { 
					heatCalculators.push_back(new MorHorn());
				}else if (globaldata->Estimators[i] == "braycurtis") { 
					heatCalculators.push_back(new BrayCurtis());
				}
			}
		}
		
		//reset calc for next command
		globaldata->setCalc("");


	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the HeatMapSimCommand class Function HeatMapSimCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the HeatMapSimCommand class function HeatMapSimCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}
//**********************************************************************************************************************

HeatMapSimCommand::~HeatMapSimCommand(){
	delete input;
	delete read;
	delete heatmap;
}

//**********************************************************************************************************************

int HeatMapSimCommand::execute(){
	try {
		int count = 1;	
		
		//if the users entered no valid calculators don't execute command
		if (heatCalculators.size() == 0) { cout << "No valid calculators." << endl; return 0; }
		
		//you have groups
		read = new ReadOTUFile(globaldata->inputFileName);	
		read->read(&*globaldata); 
			
		input = globaldata->ginput;
		lookup = input->getSharedRAbundVectors();
		vector<SharedRAbundVector*> lastLookup = lookup;
		
		if (lookup.size() < 2) { cout << "You have not provided enough valid groups.  I cannot run the command." << endl; return 0;}
				
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = globaldata->labels;
		set<int> userLines = globaldata->lines;

		
		//as long as you are not at the end of the file or done wih the lines you want
		while((lookup[0] != NULL) && ((globaldata->allLines == 1) || (userLabels.size() != 0) || (userLines.size() != 0))) {
		
			if(globaldata->allLines == 1 || globaldata->lines.count(count) == 1 || globaldata->labels.count(lookup[0]->getLabel()) == 1){			
	
				cout << lookup[0]->getLabel() << '\t' << count << endl;
				heatmap->getPic(lookup, heatCalculators);
					
				processedLabels.insert(lookup[0]->getLabel());
				userLabels.erase(lookup[0]->getLabel());
				userLines.erase(count);
			}
				
			if ((anyLabelsToProcess(lookup[0]->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLookup[0]->getLabel()) != 1)) {
				cout << lastLookup[0]->getLabel() << '\t' << count << endl;
				heatmap->getPic(lastLookup, heatCalculators);
					
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
			heatmap->getPic(lastLookup, heatCalculators);
		}
		
		for (int i = 0; i < lastLookup.size(); i++) {  delete lastLookup[i];  }
			
		//reset groups parameter
		globaldata->Groups.clear();  
		globaldata->setGroups("");
		
		return 0;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the HeatMapSimCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the HeatMapSimCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		
}

//**********************************************************************************************************************
