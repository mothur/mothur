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

GetRAbundCommand::GetRAbundCommand(string option)  {
	try {
		globaldata = GlobalData::getInstance();
		abort = false;
		allLines = 1;
		labels.clear();
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"label","sorted","outputdir","inputdir"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			
			//check to make sure all parameters are valid for command
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			string outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += hasPath(globaldata->inputFileName); //if user entered a file with a path then preserve it	
			}
			
			//make sure the user has already run the read.otu command
			if (globaldata->getListFile() == "") { m->mothurOut("You must read a listfile before you can use the get.rabund command."); m->mothurOutEndLine(); abort = true; }
			
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
				filename = outputDir + getRootName(getSimpleName(globaldata->inputFileName)) + "rabund";
				openOutputFile(filename, out);
			}
		}

	}
	catch(exception& e) {
		m->errorOut(e, "GetRAbundCommand", "GetRAbundCommand");
		exit(1);
	}			
}
//**********************************************************************************************************************

void GetRAbundCommand::help(){
	try {
		m->mothurOut("The get.rabund command can only be executed after a successful read.otu of a listfile.\n");
		m->mothurOut("The get.rabund command parameters are label and sorted.  No parameters are required.\n");
		m->mothurOut("The label parameter allows you to select what distance levels you would like included in your .rabund file, and are separated by dashes.\n");
		m->mothurOut("The sorted parameters allows you to print the rabund results sorted by abundance or not.  The default is sorted.\n");
		m->mothurOut("The get.rabund command should be in the following format: get.rabund(label=yourLabels, sorted=yourSorted).\n");
		m->mothurOut("Example get.rabund(sorted=F).\n");
		m->mothurOut("The default value for label is all labels in your inputfile.\n");
		m->mothurOut("The get.rabund command outputs a .rabund file containing the lines you selected.\n");
		m->mothurOut("Note: No spaces between parameter labels (i.e. label), '=' and parameters (i.e.yourLabels).\n\n");
	}
	catch(exception& e) {
		m->errorOut(e, "GetRAbundCommand", "help");
		exit(1);
	}
}

//**********************************************************************************************************************

GetRAbundCommand::~GetRAbundCommand(){}

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
		
		if (m->control_pressed) {  out.close(); remove(filename.c_str());  delete list; globaldata->gListVector = NULL;  return 0; }
		
		while((list != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			
			if(allLines == 1 || labels.count(list->getLabel()) == 1){
					m->mothurOut(list->getLabel()); m->mothurOutEndLine();
					rabund = new RAbundVector();				
					*rabund = (list->getRAbundVector());
					
					if (m->control_pressed) {  out.close(); remove(filename.c_str());  delete list; delete rabund; globaldata->gListVector = NULL;  return 0; }

					
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
					
					m->mothurOut(list->getLabel()); m->mothurOutEndLine();
					rabund = new RAbundVector();
					*rabund = (list->getRAbundVector());
					
					if (m->control_pressed) {  out.close(); remove(filename.c_str());  delete list; delete rabund; globaldata->gListVector = NULL;  return 0; }
					
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
			m->mothurOut("Your file does not include the label " + *it); 
			if (processedLabels.count(lastLabel) != 1) {
				m->mothurOut(". I will use " + lastLabel + "."); m->mothurOutEndLine();
				needToRun = true;
			}else {
				m->mothurOut(". Please refer to " + lastLabel + "."); m->mothurOutEndLine();
			}
		}
		
		//run last label if you need to
		if (needToRun == true)  {
			if (list != NULL) {	delete list;	}
			list = input->getListVector(lastLabel);
			
			m->mothurOut(list->getLabel()); m->mothurOutEndLine();
			rabund = new RAbundVector();
			*rabund = (list->getRAbundVector());
			
			if (m->control_pressed) {  out.close(); remove(filename.c_str());  delete list; delete rabund; globaldata->gListVector = NULL;  return 0; }
			
			if(sorted)	{   rabund->print(out);				}
			else		{	rabund->nonSortedPrint(out);	}

			delete rabund;
			delete list;
		}
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Name: "); m->mothurOutEndLine();
		m->mothurOut(filename); m->mothurOutEndLine();	
		m->mothurOutEndLine();
		
		out.close(); 
		
		globaldata->gListVector = NULL;
		
		return 0;		
	}

	catch(exception& e) {
		m->errorOut(e, "GetRAbundCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************


