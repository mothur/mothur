/*
 *  getsabundcommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 6/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "getsabundcommand.h"

//**********************************************************************************************************************

GetSAbundCommand::GetSAbundCommand(){
	try {
		globaldata = GlobalData::getInstance();
		filename = getRootName(globaldata->inputFileName) + "sabund";
		
		openOutputFile(filename, out);
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the GetSAbundCommand class Function GetSAbundCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the GetSAbundCommand class function GetSAbundCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
			
}

//**********************************************************************************************************************

GetSAbundCommand::~GetSAbundCommand(){
}

//**********************************************************************************************************************

int GetSAbundCommand::execute(){
	try {
		int count = 1;
		
		//using order vector so you don't have to distinguish between the list and rabund files
		read = new ReadOTUFile(globaldata->inputFileName);	
		read->read(&*globaldata); 
		
		order = globaldata->gorder;
		lastOrder = order;
		input = globaldata->ginput;
						
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = globaldata->labels;
		set<int> userLines = globaldata->lines;

		
		while((order != NULL) && ((globaldata->allLines == 1) || (userLabels.size() != 0) || (userLines.size() != 0))) {
			
			if(globaldata->allLines == 1 || globaldata->lines.count(count) == 1 || globaldata->labels.count(order->getLabel()) == 1){
					cout << order->getLabel() << '\t' << count << endl;
					sabund = new SAbundVector();
					*sabund = (order->getSAbundVector());
					sabund->print(out);
					delete sabund;

					processedLabels.insert(order->getLabel());
					userLabels.erase(order->getLabel());
					userLines.erase(count);
			}
			
			if ((anyLabelsToProcess(order->getLabel(), userLabels, "") == true) && (processedLabels.count(lastOrder->getLabel()) != 1)) {
					cout << lastOrder->getLabel() << '\t' << count << endl;
					sabund = new SAbundVector();
					*sabund = (lastOrder->getSAbundVector());
					sabund->print(out);
					delete sabund;

					processedLabels.insert(lastOrder->getLabel());
					userLabels.erase(lastOrder->getLabel());
			}
			
			if (count != 1) { delete lastOrder; }
			lastOrder = order;	
					
			order = (input->getOrderVector());
			count++;
		}
		
		//output error messages about any remaining user labels
		set<string>::iterator it;
		bool needToRun = false;
		for (it = userLabels.begin(); it != userLabels.end(); it++) {  
			cout << "Your file does not include the label "<< *it; 
			if (processedLabels.count(lastOrder->getLabel()) != 1) {
				cout << ". I will use " << lastOrder->getLabel() << "." << endl;
				needToRun = true;
			}else {
				cout << ". Please refer to " << lastOrder->getLabel() << "." << endl;
			}
		}
		
		//run last line if you need to
		if (needToRun == true)  {
			cout << lastOrder->getLabel() << '\t' << count << endl;
			sabund = new SAbundVector();
			*sabund = (lastOrder->getSAbundVector());
			sabund->print(out);
			delete sabund;
		}
		delete lastOrder;

		out.close();
		return 0;		
	}

	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the GetSAbundCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the GetSAbundCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

//**********************************************************************************************************************


