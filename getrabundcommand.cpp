/*
 *  getrabundcommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 6/2/09.
 *  Copyright 2009 Schloss Lab Umass Amherst. All rights reserved.
 *
 */

#include "getrabundcommand.h"

//**********************************************************************************************************************

GetRAbundCommand::GetRAbundCommand(string option){
	try {
		globaldata = GlobalData::getInstance();
		abort = false;
		allLines = 1;
		labels.clear();
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"label","sorted"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			
			//check to make sure all parameters are valid for command
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//make sure the user has already run the read.otu command
			if (globaldata->getListFile() == "") { mothurOut("You must read a listfile before you can use the get.rabund command."); mothurOutEndLine(); abort = true; }
			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			
			string temp;
			temp = validParameter.validFile(parameters, "sorted", false);			if (temp == "not found") { temp = "T"; }
			sorted = isTrue(temp);
			
			label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
			
			//if the user has not specified any labels use the ones from read.otu
			if(label == "") {  
				allLines = globaldata->allLines; 
				labels = globaldata->labels; 
			}
				
			if (abort == false) {
				filename = getRootName(globaldata->inputFileName) + "rabund";
				openOutputFile(filename, out);
			}
		}

	}
	catch(exception& e) {
		errorOut(e, "GetRAbundCommand", "GetRAbundCommand");
		exit(1);
	}			
}
//**********************************************************************************************************************

void GetRAbundCommand::help(){
	try {
		mothurOut("The get.rabund command can only be executed after a successful read.otu of a listfile.\n");
		mothurOut("The get.rabund command parameters are label and sorted.  No parameters are required.\n");
		mothurOut("The label parameter allows you to select what distance levels you would like included in your .rabund file, and are separated by dashes.\n");
		mothurOut("The sorted parameters allows you to print the rabund results sorted by abundance or not.  The default is sorted.\n");
		mothurOut("The get.rabund command should be in the following format: get.rabund(label=yourLabels, sorted=yourSorted).\n");
		mothurOut("Example get.rabund(sorted=F).\n");
		mothurOut("The default value for label is all labels in your inputfile.\n");
		mothurOut("The get.rabund command outputs a .rabund file containing the lines you selected.\n");
		mothurOut("Note: No spaces between parameter labels (i.e. label), '=' and parameters (i.e.yourLabels).\n\n");
	}
	catch(exception& e) {
		errorOut(e, "GetRAbundCommand", "help");
		exit(1);
	}
}

//**********************************************************************************************************************

GetRAbundCommand::~GetRAbundCommand(){
	if (abort == false) {  globaldata->gListVector = NULL; }
}

//**********************************************************************************************************************

int GetRAbundCommand::execute(){
	try {
	
		if (abort == true) { return 0; }
		
		//read first line
		read = new ReadOTUFile(globaldata->inputFileName);	
		read->read(&*globaldata); 
			
		input = globaldata->ginput;
		list = globaldata->gListVector;
		string lastLabel = list->getLabel();
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;
		
		while((list != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			
			if(allLines == 1 || labels.count(list->getLabel()) == 1){
					mothurOut(list->getLabel()); mothurOutEndLine();
					rabund = new RAbundVector();				
					*rabund = (list->getRAbundVector());
					
					if(sorted)	{   rabund->print(out);				}
					else		{	rabund->nonSortedPrint(out);	}
					
					delete rabund;
															
					processedLabels.insert(list->getLabel());
					userLabels.erase(list->getLabel());
			}
			
			if ((anyLabelsToProcess(list->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
					string saveLabel = list->getLabel();
					
					delete list;
					list = input->getListVector(lastLabel);
					
					mothurOut(list->getLabel()); mothurOutEndLine();
					rabund = new RAbundVector();
					*rabund = (list->getRAbundVector());
					
					if(sorted)	{   rabund->print(out);				}
					else		{	rabund->nonSortedPrint(out);	}

					delete rabund;

					processedLabels.insert(list->getLabel());
					userLabels.erase(list->getLabel());
					
					//restore real lastlabel to save below
					list->setLabel(saveLabel);
			}
			
			lastLabel = list->getLabel();		
			
			delete list;
			list = input->getListVector();
		}
		
		//output error messages about any remaining user labels
		set<string>::iterator it;
		bool needToRun = false;
		for (it = userLabels.begin(); it != userLabels.end(); it++) {  
			mothurOut("Your file does not include the label " + *it); 
			if (processedLabels.count(lastLabel) != 1) {
				mothurOut(". I will use " + lastLabel + "."); mothurOutEndLine();
				needToRun = true;
			}else {
				mothurOut(". Please refer to " + lastLabel + "."); mothurOutEndLine();
			}
		}
		
		//run last label if you need to
		if (needToRun == true)  {
			if (list != NULL) {	delete list;	}
			list = input->getListVector(lastLabel);
			
			mothurOut(list->getLabel()); mothurOutEndLine();
			rabund = new RAbundVector();
			*rabund = (list->getRAbundVector());
			
			if(sorted)	{   rabund->print(out);				}
			else		{	rabund->nonSortedPrint(out);	}

			delete rabund;
			delete list;
		}

		out.close(); 
		return 0;		
	}

	catch(exception& e) {
		errorOut(e, "GetRAbundCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************


