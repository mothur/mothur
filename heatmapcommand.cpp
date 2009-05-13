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
	
		if (format == "sharedfile") {
			//you have groups
			read = new ReadOTUFile(globaldata->inputFileName);	
			read->read(&*globaldata); 
			
			input = globaldata->ginput;
			lookup = input->getSharedRAbundVectors();
		}else if (format == "list") {
			//you are using just a list file and have only one group
			read = new ReadOTUFile(globaldata->inputFileName);	
			read->read(&*globaldata); 
			
			rabund = globaldata->rabund;
			input = globaldata->ginput;		
		}
		
		if (format != "list") {	
		
			while(lookup[0] != NULL){
		
				if(globaldata->allLines == 1 || globaldata->lines.count(count) == 1 || globaldata->labels.count(lookup[0]->getLabel()) == 1){			
	
					cout << lookup[0]->getLabel() << '\t' << count << endl;
					heatmap->getPic(lookup);
				}
						
				//get next line to process
				lookup = input->getSharedRAbundVectors();				
				count++;
			}
			
			//reset groups parameter
			globaldata->Groups.clear();  
			
		}else{
		
			while(rabund != NULL){
		
				if(globaldata->allLines == 1 || globaldata->lines.count(count) == 1 || globaldata->labels.count(rabund->getLabel()) == 1){			
	
					cout << rabund->getLabel() << '\t' << count << endl;
					heatmap->getPic(rabund);
				}
				
				rabund = input->getRAbundVector();
				count++;
			}
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


