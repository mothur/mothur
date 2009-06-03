/*
 *  heatmapcommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 3/25/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "heatmapcommand.h"


//**********************************************************************************************************************

HeatMapCommand::HeatMapCommand(){
	try {
		globaldata = GlobalData::getInstance();
		heatmap = new HeatMap();
		format = globaldata->getFormat();
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the HeatMapCommand class Function HeatMapCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the HeatMapCommand class function HeatMapCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}
//**********************************************************************************************************************

HeatMapCommand::~HeatMapCommand(){
	delete input;
	delete read;
	delete heatmap;
}

//**********************************************************************************************************************

int HeatMapCommand::execute(){
	try {
		int count = 1;	
		RAbundVector* lastRAbund;
		vector<SharedRAbundVector*> lastLookup;
	
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
			
			rabund = globaldata->rabund;
			lastRAbund = globaldata->rabund;
			input = globaldata->ginput;		
		}
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = globaldata->labels;
		set<int> userLines = globaldata->lines;

		if (format != "list") {	
		
			//as long as you are not at the end of the file or done wih the lines you want
			while((lookup[0] != NULL) && ((globaldata->allLines == 1) || (userLabels.size() != 0) || (userLines.size() != 0))) {
		
				if(globaldata->allLines == 1 || globaldata->lines.count(count) == 1 || globaldata->labels.count(lookup[0]->getLabel()) == 1){			
	
					cout << lookup[0]->getLabel() << '\t' << count << endl;
					heatmap->getPic(lookup);
					
					processedLabels.insert(lookup[0]->getLabel());
					userLabels.erase(lookup[0]->getLabel());
					userLines.erase(count);
				}
				
				if ((anyLabelsToProcess(lookup[0]->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLookup[0]->getLabel()) != 1)) {
					cout << lastLookup[0]->getLabel() << '\t' << count << endl;
					heatmap->getPic(lastLookup);
					
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
				heatmap->getPic(lastLookup);
			}
		
			for (int i = 0; i < lastLookup.size(); i++) {  delete lastLookup[i];  }
			
			//reset groups parameter
			globaldata->Groups.clear();  
			
		}else{
		
			while((rabund != NULL) && ((globaldata->allLines == 1) || (userLabels.size() != 0) || (userLines.size() != 0))) {

				if(globaldata->allLines == 1 || globaldata->lines.count(count) == 1 || globaldata->labels.count(rabund->getLabel()) == 1){			
	
					cout << rabund->getLabel() << '\t' << count << endl;
					heatmap->getPic(rabund);
					
					processedLabels.insert(rabund->getLabel());
					userLabels.erase(rabund->getLabel());
					userLines.erase(count);
				}
				
				if ((anyLabelsToProcess(rabund->getLabel(), userLabels, "") == true) && (processedLabels.count(lastRAbund->getLabel()) != 1)) {

					cout << lastRAbund->getLabel() << '\t' << count << endl;
					heatmap->getPic(lastRAbund);
					
					processedLabels.insert(lastRAbund->getLabel());
					userLabels.erase(lastRAbund->getLabel());
				}		
				
				if (count != 1) { delete lastRAbund; }
				lastRAbund = rabund;			

				rabund = input->getRAbundVector();
				count++;
			}
			
			//output error messages about any remaining user labels
			set<string>::iterator it;
			bool needToRun = false;
			for (it = userLabels.begin(); it != userLabels.end(); it++) {  
				cout << "Your file does not include the label "<< *it; 
				if (processedLabels.count(lastRAbund->getLabel()) != 1) {
					cout << ". I will use " << lastRAbund->getLabel() << "." << endl;
					needToRun = true;
				}else {
					cout << ". Please refer to " << lastRAbund->getLabel() << "." << endl;
				}
			}
		
			//run last line if you need to
			if (needToRun == true)  {
				cout << lastRAbund->getLabel() << '\t' << count << endl;
				heatmap->getPic(lastRAbund);
			}
		
			delete lastRAbund;

		}
		
		globaldata->setGroups("");
		return 0;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the HeatMapCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the HeatMapCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		
}

//**********************************************************************************************************************


