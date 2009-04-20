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
		util = new SharedUtil();
		
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
	delete util;
}

//**********************************************************************************************************************

int HeatMapCommand::execute(){
	try {
		int count = 1;	
		
		if (format == "sharedfile") {
			//you have groups
			read = new ReadPhilFile(globaldata->inputFileName);	
			read->read(&*globaldata); 
			
			input = globaldata->ginput;
			order = input->getSharedOrderVector();
		}else if (format == "shared") {
			//you are using a list and a groupfile
			read = new ReadPhilFile(globaldata->inputFileName);	
			read->read(&*globaldata); 
		
			input = globaldata->ginput;
			SharedList = globaldata->gSharedList;
			order = SharedList->getSharedOrderVector();
		}else if (format == "list") {
			//you are using just a list file and have only one group
			read = new ReadPhilFile(globaldata->inputFileName);	
			read->read(&*globaldata); 
			
			ordersingle = globaldata->gorder;
			input = globaldata->ginput;
		}
		
		if (format != "list") {	
		
			util->setGroups(globaldata->Groups, globaldata->gGroupmap->namesOfGroups, "heat");
			globaldata->setGroups("");


			while(order != NULL){
		
				if(globaldata->allLines == 1 || globaldata->lines.count(count) == 1 || globaldata->labels.count(order->getLabel()) == 1){			
	
					cout << order->getLabel() << '\t' << count << endl;
					heatmap->getPic(order);

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
			globaldata->Groups.clear();  
			
		}else{
			while(ordersingle != NULL){
		
				if(globaldata->allLines == 1 || globaldata->lines.count(count) == 1 || globaldata->labels.count(ordersingle->getLabel()) == 1){			
	
					cout << ordersingle->getLabel() << '\t' << count << endl;
					heatmap->getPic(ordersingle);
					
				}
				
				ordersingle = (input->getOrderVector());
				count++;
			}
		}
		
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


