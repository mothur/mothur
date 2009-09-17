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
		
		groupMap = globaldata->gGroupmap;
		
		//fill filehandles with neccessary ofstreams
		int i;
		ofstream* temp;
		for (i=0; i<groupMap->getNumGroups(); i++) {
			temp = new ofstream;
			filehandles[groupMap->namesOfGroups[i]] = temp;
		}
		
		//set fileroot
		fileroot = getRootName(globaldata->getListFile());
		
		//clears file before we start to write to it below
		for (int i=0; i<groupMap->getNumGroups(); i++) {
			remove((fileroot + groupMap->namesOfGroups[i] + ".rabund").c_str());
		}

	}
	catch(exception& e) {
		errorOut(e, "SharedCommand", "SharedCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int SharedCommand::execute(){
	try {
		
		//lookup.clear();
		int count = 1;
		string errorOff = "no error";
			
		//read in listfile
		read = new ReadOTUFile(globaldata->inputFileName);	
		read->read(&*globaldata); 
		delete read;

		input = globaldata->ginput;
		SharedList = globaldata->gSharedList;
		string lastLabel = SharedList->getLabel();
		vector<SharedRAbundVector*> lookup; 
		
		if (SharedList->getNumSeqs() != groupMap->getNumSeqs()) {  
			mothurOut("Your group file contains " + toString(groupMap->getNumSeqs()) + " sequences and list file contains " + toString(SharedList->getNumSeqs()) + " sequences. Please correct."); mothurOutEndLine(); 
			
			//delete memory
			for (it3 = filehandles.begin(); it3 != filehandles.end(); it3++) {
				delete it3->second;
			}
			delete SharedList;
			globaldata->gSharedList = NULL;
			
			return(0); 
		}
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = globaldata->labels;
		set<int> userLines = globaldata->lines;
		
		
		while((SharedList != NULL) && ((globaldata->allLines == 1) || (userLabels.size() != 0) || (userLines.size() != 0))) {
			

			if(globaldata->allLines == 1 || globaldata->lines.count(count) == 1 || globaldata->labels.count(SharedList->getLabel()) == 1){
					
					lookup = SharedList->getSharedRAbundVector();
					mothurOut(lookup[0]->getLabel()); mothurOutEndLine();
					
					printSharedData(lookup); //prints info to the .shared file
					for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }
				
					processedLabels.insert(SharedList->getLabel());
					userLabels.erase(SharedList->getLabel());
					userLines.erase(count);
			}
			
			if ((anyLabelsToProcess(SharedList->getLabel(), userLabels, errorOff) == true) && (processedLabels.count(lastLabel) != 1)) {
					delete SharedList;
					SharedList = input->getSharedListVector(lastLabel); //get new list vector to process
					
					lookup = SharedList->getSharedRAbundVector();
					mothurOut(lookup[0]->getLabel()); mothurOutEndLine();
					
					printSharedData(lookup); //prints info to the .shared file
					for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }
					
					processedLabels.insert(SharedList->getLabel());
					userLabels.erase(SharedList->getLabel());
			}
			
		
			lastLabel = SharedList->getLabel();
				
			delete SharedList;
			SharedList = input->getSharedListVector(); //get new list vector to process
			
			count++;		
		}
		
		//output error messages about any remaining user labels
		set<string>::iterator it;
		bool needToRun = false;
		for (it = userLabels.begin(); it != userLabels.end(); it++) {  
			if (processedLabels.count(lastLabel) != 1) {
				needToRun = true;
			}
		}
		
		//run last line if you need to
		if (needToRun == true)  {
			if (SharedList != NULL) {	delete SharedList;	}
			SharedList = input->getSharedListVector(lastLabel); //get new list vector to process
					
			lookup = SharedList->getSharedRAbundVector();
			mothurOut(lookup[0]->getLabel()); mothurOutEndLine();
			
			printSharedData(lookup); //prints info to the .shared file
			for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }
			delete SharedList;
		}
		
		globaldata->gSharedList = NULL;
		
		out.close();
		
		for (it3 = filehandles.begin(); it3 != filehandles.end(); it3++) {
			delete it3->second;
		}

		
		return 0;
	}
	catch(exception& e) {
		errorOut(e, "SharedCommand", "execute");
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
			
			RAbundVector rav = thislookup[i]->getRAbundVector();
			openOutputFileAppend(fileroot + thislookup[i]->getGroup() + ".rabund", *(filehandles[thislookup[i]->getGroup()]));
			rav.print(*(filehandles[thislookup[i]->getGroup()]));
			(*(filehandles[thislookup[i]->getGroup()])).close();
		}
 
	}
	catch(exception& e) {
		errorOut(e, "SharedCommand", "printSharedData");
		exit(1);
	}
}

//**********************************************************************************************************************

SharedCommand::~SharedCommand(){
	//delete list;
	
	
}

//**********************************************************************************************************************
