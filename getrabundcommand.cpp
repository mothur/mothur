/*
 *  getrabundcommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 6/2/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "getrabundcommand.h"

//**********************************************************************************************************************

GetRAbundCommand::GetRAbundCommand(){
	try {
		globaldata = GlobalData::getInstance();
		filename = getRootName(globaldata->inputFileName) + "rabund";
		
		openOutputFile(filename, out);
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the GetRAbundCommand class Function GetRAbundCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the GetRAbundCommand class function GetRAbundCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
			
}

//**********************************************************************************************************************

GetRAbundCommand::~GetRAbundCommand(){
}

//**********************************************************************************************************************

int GetRAbundCommand::execute(){
	try {
		int count = 1;
		
		//read first line
		read = new ReadOTUFile(globaldata->inputFileName);	
		read->read(&*globaldata); 
			
		input = globaldata->ginput;
		list = globaldata->gListVector;
		ListVector* lastList = list;
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = globaldata->labels;
		set<int> userLines = globaldata->lines;

		
		while((list != NULL) && ((globaldata->allLines == 1) || (userLabels.size() != 0) || (userLines.size() != 0))) {
			
			if(globaldata->allLines == 1 || globaldata->lines.count(count) == 1 || globaldata->labels.count(list->getLabel()) == 1){
					cout << list->getLabel() << '\t' << count << endl;
					rabund = new RAbundVector();
					*rabund = (list->getRAbundVector());
					rabund->print(out);
					delete rabund;
															
					processedLabels.insert(list->getLabel());
					userLabels.erase(list->getLabel());
					userLines.erase(count);
			}
			
			if ((anyLabelsToProcess(list->getLabel(), userLabels, "") == true) && (processedLabels.count(lastList->getLabel()) != 1)) {
					cout << lastList->getLabel() << '\t' << count << endl;
					rabund = new RAbundVector();
					*rabund = (lastList->getRAbundVector());
					rabund->print(out);
					delete rabund;

					processedLabels.insert(lastList->getLabel());
					userLabels.erase(lastList->getLabel());
			}
			
			if (count != 1) { delete lastList; }
			lastList = list;			
			
			list = input->getListVector();
			count++;
		}
		
		//output error messages about any remaining user labels
		set<string>::iterator it;
		bool needToRun = false;
		for (it = userLabels.begin(); it != userLabels.end(); it++) {  
			cout << "Your file does not include the label "<< *it; 
			if (processedLabels.count(lastList->getLabel()) != 1) {
				cout << ". I will use " << lastList->getLabel() << "." << endl;
				needToRun = true;
			}else {
				cout << ". Please refer to " << lastList->getLabel() << "." << endl;
			}
		}
		
		//run last line if you need to
		if (needToRun == true)  {
			cout << lastList->getLabel() << '\t' << count << endl;
			rabund = new RAbundVector();
			*rabund = (lastList->getRAbundVector());
			rabund->print(out);
			delete rabund;
		}
		delete lastList;

		out.close(); 
		return 0;		
	}

	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the GetRAbundCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the GetRAbundCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

//**********************************************************************************************************************


