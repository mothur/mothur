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
vector<string> GetRAbundCommand::setParameters(){	
	try {
		CommandParameter plist("list", "InputTypes", "", "", "LRSS", "LRSS", "none",false,false); parameters.push_back(plist);
		CommandParameter psabund("sabund", "InputTypes", "", "", "LRSS", "LRSS", "none",false,false); parameters.push_back(psabund);		
		CommandParameter psorted("sorted", "Boolean", "", "T", "", "", "",false,false); parameters.push_back(psorted);
		CommandParameter plabel("label", "String", "", "", "", "", "",false,false); parameters.push_back(plabel);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "GetRAbundCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string GetRAbundCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The get.rabund command parameters are list, sabund, label and sorted.  list or sabund parameters are required, unless you have valid current files.\n";
		helpString += "The label parameter allows you to select what distance levels you would like included in your .rabund file, and are separated by dashes.\n";
		helpString += "The sorted parameters allows you to print the rabund results sorted by abundance or not.  The default is sorted.\n";
		helpString += "The get.rabund command should be in the following format: get.rabund(label=yourLabels, sorted=yourSorted).\n";
		helpString += "Example get.rabund(sorted=F).\n";
		helpString += "The default value for label is all labels in your inputfile.\n";
		helpString += "The get.rabund command outputs a .rabund file containing the lines you selected.\n";
		helpString += "Note: No spaces between parameter labels (i.e. label), '=' and parameters (i.e.yourLabels).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "GetRAbundCommand", "getHelpString");
		exit(1);
	}
}

//**********************************************************************************************************************
GetRAbundCommand::GetRAbundCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["rabund"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "GetRAbundCommand", "GetRAbundCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
GetRAbundCommand::GetRAbundCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		allLines = 1;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			map<string,string>::iterator it;
			
			ValidParameters validParameter;
			
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["rabund"] = tempOutNames;
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("list");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["list"] = inputDir + it->second;		}
				}
				
				it = parameters.find("sabund");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["sabund"] = inputDir + it->second;		}
				}
			}
			
			
			//check for required parameters
			listfile = validParameter.validFile(parameters, "list", true);
			if (listfile == "not open") { listfile = ""; abort = true; }
			else if (listfile == "not found") { listfile = ""; }
			else {  format = "list"; inputfile = listfile; m->setListFile(listfile); }
			
			sabundfile = validParameter.validFile(parameters, "sabund", true);
			if (sabundfile == "not open") { sabundfile = ""; abort = true; }	
			else if (sabundfile == "not found") { sabundfile = ""; }
			else {  format = "sabund"; inputfile = sabundfile; m->setSabundFile(sabundfile); }
			
			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			string temp;
			temp = validParameter.validFile(parameters, "sorted", false);			if (temp == "not found") { temp = "T"; }
			sorted = m->isTrue(temp);
			
			label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  m->splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
			
			if ((listfile == "") && (sabundfile == "")) { 
				//is there are current file available for any of these?
				//give priority to shared, then list, then rabund, then sabund
				//if there is a current shared file, use it
				listfile = m->getListFile(); 
				if (listfile != "") { inputfile = listfile; format = "list"; m->mothurOut("Using " + listfile + " as input file for the list parameter."); m->mothurOutEndLine(); }
				else { 
					sabundfile = m->getSabundFile(); 
					if (sabundfile != "") { inputfile = sabundfile; format = "sabund"; m->mothurOut("Using " + sabundfile + " as input file for the sabund parameter."); m->mothurOutEndLine(); }
					else { 
						m->mothurOut("No valid current files. You must provide a list or sabund file."); m->mothurOutEndLine(); 
						abort = true;
					}
				}
			}
			
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = m->hasPath(inputfile); 	}			
			
		}
			

	}
	catch(exception& e) {
		m->errorOut(e, "GetRAbundCommand", "GetRAbundCommand");
		exit(1);
	}			
}
//**********************************************************************************************************************

int GetRAbundCommand::execute(){
	try {
	
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		filename = outputDir + m->getRootName(m->getSimpleName(inputfile)) + "rabund";
		m->openOutputFile(filename, out);
		
		input = new InputData(inputfile, format);
		rabund = input->getRAbundVector();
		string lastLabel = rabund->getLabel();
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;
		
		if (m->control_pressed) {  outputTypes.clear();  out.close(); m->mothurRemove(filename); delete rabund; delete input; return 0; }
		
		while((rabund != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			
			if(allLines == 1 || labels.count(rabund->getLabel()) == 1){
					m->mothurOut(rabund->getLabel()); m->mothurOutEndLine();
					
					if (m->control_pressed) {   outputTypes.clear(); out.close(); m->mothurRemove(filename);  delete input; delete rabund;  return 0;  }
					
					if(sorted)	{   rabund->print(out);				}
					else		{	rabund->nonSortedPrint(out);	}
															
					processedLabels.insert(rabund->getLabel());
					userLabels.erase(rabund->getLabel());
			}
			
			if ((m->anyLabelsToProcess(rabund->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
					string saveLabel = rabund->getLabel();
					
					delete rabund;
					rabund = input->getRAbundVector(lastLabel);
					
					m->mothurOut(rabund->getLabel()); m->mothurOutEndLine();
					
					if (m->control_pressed) {   outputTypes.clear(); out.close(); m->mothurRemove(filename);  delete input; delete rabund;  return 0;  }
					
					if(sorted)	{   rabund->print(out);				}
					else		{	rabund->nonSortedPrint(out);	}

					processedLabels.insert(rabund->getLabel());
					userLabels.erase(rabund->getLabel());
					
					//restore real lastlabel to save below
					rabund->setLabel(saveLabel);
			}
			
			lastLabel = rabund->getLabel();		
			
			delete rabund;
			rabund = input->getRAbundVector();
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
			if (rabund != NULL) {	delete rabund;	}
			rabund = input->getRAbundVector(lastLabel);
			
			m->mothurOut(rabund->getLabel()); m->mothurOutEndLine();
					
			if (m->control_pressed) {  outputTypes.clear(); out.close(); m->mothurRemove(filename);  delete input; delete rabund;  return 0; }
			
			if(sorted)	{   rabund->print(out);				}
			else		{	rabund->nonSortedPrint(out);	}

			delete rabund;
		}
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Name: "); m->mothurOutEndLine();
		m->mothurOut(filename); m->mothurOutEndLine();	outputNames.push_back(filename); outputTypes["rabund"].push_back(filename);
		m->mothurOutEndLine();
		
		out.close(); 
				
		//set rabund file as new current rabundfile
		string current = "";
		itTypes = outputTypes.find("rabund");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setRabundFile(current); }
		}
		
		return 0;		
	}

	catch(exception& e) {
		m->errorOut(e, "GetRAbundCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************


