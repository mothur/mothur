/*
 *  sharedcommand.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "sharedcommand.h"

//**********************************************************************************************************************

SharedCommand::SharedCommand(){
	try {
		globaldata = GlobalData::getInstance();
		
		//getting output filename
		filename = globaldata->inputFileName;
		filename = getRootName(filename);
		filename = filename + "shared";
		openOutputFile(filename, out);
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SharedCommand class Function SharedCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SharedCommand class function SharedCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}

}
//**********************************************************************************************************************

int SharedCommand::execute(){
	try {
		
		cout << "creating sharedfile...";
		//lookup.clear();
		int count = 1;
		string errorOff = "no error";
			
		//read in listfile
		read = new ReadOTUFile(globaldata->inputFileName);	
		read->read(&*globaldata); 

		input = globaldata->ginput;
		SharedList = globaldata->gSharedList;
		SharedListVector* lastList = SharedList;
		vector<SharedRAbundVector*> lookup; 
				
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = globaldata->labels;
		set<int> userLines = globaldata->lines;
		
		
		while((SharedList != NULL) && ((globaldata->allLines == 1) || (userLabels.size() != 0) || (userLines.size() != 0))) {
			

			if(globaldata->allLines == 1 || globaldata->lines.count(count) == 1 || globaldata->labels.count(SharedList->getLabel()) == 1){
					lookup = SharedList->getSharedRAbundVector();
					printSharedData(lookup); //prints info to the .shared file
					for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }
				
					processedLabels.insert(SharedList->getLabel());
					userLabels.erase(SharedList->getLabel());
					userLines.erase(count);
			}
			
			if ((anyLabelsToProcess(SharedList->getLabel(), userLabels, errorOff) == true) && (processedLabels.count(lastList->getLabel()) != 1)) {
					lookup = lastList->getSharedRAbundVector();
					printSharedData(lookup); //prints info to the .shared file
					for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }
					
					processedLabels.insert(lastList->getLabel());
					userLabels.erase(lastList->getLabel());
			}
			
			if (count != 1) { delete lastList; }
			lastList = SharedList;	
						 		
			SharedList = input->getSharedListVector(); //get new list vector to process
			
			count++;		
		}
		
		//output error messages about any remaining user labels
		set<string>::iterator it;
		bool needToRun = false;
		for (it = userLabels.begin(); it != userLabels.end(); it++) {  
			//cout << "Your file does not include the label "<< *it; 
			if (processedLabels.count(lastList->getLabel()) != 1) {
				//cout << ". I will use " << lastList->getLabel() << "." << endl;
				needToRun = true;
			}else {
				//cout << ". Please refer to " << lastList->getLabel() << "." << endl;
			}
		}
		
		//run last line if you need to
		if (needToRun == true)  {
			lookup = lastList->getSharedRAbundVector();
			printSharedData(lookup); //prints info to the .shared file
			for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }

		}
		
		delete lastList; globaldata->gSharedList = NULL;
		delete read;
		
		out.close();
		
		cout << "complete." << endl;
		return 0;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SharedCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SharedCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}

}
//**********************************************************************************************************************
void SharedCommand::printSharedData(vector<SharedRAbundVector*> thislookup) {
	try {
		
		//initialize bin values
		for (int i = 0; i < thislookup.size(); i++) {
			out << thislookup[i]->getLabel() << '\t' << thislookup[i]->getGroup() << '\t';
			thislookup[i]->print(out);
		}
 
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the SharedCommand class Function printSharedData. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the SharedCommand class function printSharedData. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	
}

//**********************************************************************************************************************

SharedCommand::~SharedCommand(){
	//delete list;
	
	
}

//**********************************************************************************************************************
