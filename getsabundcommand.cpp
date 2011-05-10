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
vector<string> GetSAbundCommand::setParameters(){	
	try {
		CommandParameter plist("list", "InputTypes", "", "", "LRSS", "LRSS", "none",false,false); parameters.push_back(plist);
		CommandParameter prabund("rabund", "InputTypes", "", "", "LRSS", "LRSS", "none",false,false); parameters.push_back(prabund);		
		CommandParameter plabel("label", "String", "", "", "", "", "",false,false); parameters.push_back(plabel);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "GetSAbundCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string GetSAbundCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The get.sabund command parameters is list, rabund and label.  list or rabund is required unless a valid current file exists.\n";
		helpString += "The label parameter allows you to select what distance levels you would like included in your .sabund file, and are separated by dashes.\n";
		helpString += "The get.sabund command should be in the following format: get.sabund(label=yourLabels).\n";
		helpString += "Example get.sabund().\n";
		helpString += "The default value for label is all labels in your inputfile.\n";
		helpString += "The get.sabund command outputs a .sabund file containing the labels you selected.\n";
		helpString += "Note: No spaces between parameter labels (i.e. label), '=' and parameters (i.e.yourLabel).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "GetSAbundCommand", "getHelpString");
		exit(1);
	}
}

//**********************************************************************************************************************
GetSAbundCommand::GetSAbundCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["sabund"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "GetSAbundCommand", "GetSAbundCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
GetSAbundCommand::GetSAbundCommand(string option)  {
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
			outputTypes["sabund"] = tempOutNames;
			
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
				
				it = parameters.find("rabund");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["rabund"] = inputDir + it->second;		}
				}
			}
			
			
			//check for required parameters
			listfile = validParameter.validFile(parameters, "list", true);
			if (listfile == "not open") { listfile = ""; abort = true; }
			else if (listfile == "not found") { listfile = ""; }
			else {  format = "list"; inputfile = listfile; }
			
			rabundfile = validParameter.validFile(parameters, "rabund", true);
			if (rabundfile == "not open") { rabundfile = ""; abort = true; }	
			else if (rabundfile == "not found") { rabundfile = ""; }
			else {  format = "rabund"; inputfile = rabundfile; }
			
		
						//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  m->splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
			
			if ((listfile == "") && (rabundfile == "")) { 
				//is there are current file available for any of these?
				//give priority to shared, then list, then rabund, then sabund
				//if there is a current shared file, use it
				listfile = m->getListFile(); 
				if (listfile != "") { inputfile = listfile; format = "list"; m->mothurOut("Using " + listfile + " as input file for the list parameter."); m->mothurOutEndLine(); }
				else { 
					rabundfile = m->getRabundFile(); 
					if (rabundfile != "") { inputfile = rabundfile; format = "rabund"; m->mothurOut("Using " + rabundfile + " as input file for the rabund parameter."); m->mothurOutEndLine(); }
					else { 
						m->mothurOut("No valid current files. You must provide a list or rabund file."); m->mothurOutEndLine(); 
						abort = true;
					}
				}
			}
			
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = m->hasPath(inputfile); 	}			
			
			
		}

	}
	catch(exception& e) {
		m->errorOut(e, "GetSAbundCommand", "GetSAbundCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int GetSAbundCommand::execute(){
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		filename = outputDir + m->getRootName(m->getSimpleName(inputfile)) + "sabund";
		m->openOutputFile(filename, out);
		
		input = new InputData(inputfile, format);
		sabund = input->getSAbundVector();
		string lastLabel = sabund->getLabel();
		
						
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;
		
		if (m->control_pressed) {  outputTypes.clear(); out.close(); remove(filename.c_str());  delete sabund; delete input; return 0; }

		
		while((sabund != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			
			if(allLines == 1 || labels.count(sabund->getLabel()) == 1){
					m->mothurOut(sabund->getLabel());  m->mothurOutEndLine();
					
					sabund->print(out);
					
				if (m->control_pressed) { outputTypes.clear();  out.close(); remove(filename.c_str());  delete sabund; delete input;  return 0; }

					processedLabels.insert(sabund->getLabel());
					userLabels.erase(sabund->getLabel());
			}
			
			if ((m->anyLabelsToProcess(sabund->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
					string saveLabel = sabund->getLabel();
					
					delete sabund;		
					sabund = (input->getSAbundVector(lastLabel));
					
					m->mothurOut(sabund->getLabel());  m->mothurOutEndLine();
					sabund->print(out);
					
					if (m->control_pressed) {  outputTypes.clear(); out.close(); remove(filename.c_str());  delete sabund; delete input;  return 0; }

					processedLabels.insert(sabund->getLabel());
					userLabels.erase(sabund->getLabel());
					
					//restore real lastlabel to save below
					sabund->setLabel(saveLabel);
			}
			
			
			lastLabel = sabund->getLabel();	
			
			delete sabund;		
			sabund = (input->getSAbundVector());
		}
		
		//output error messages about any remaining user labels
		set<string>::iterator it;
		bool needToRun = false;
		for (it = userLabels.begin(); it != userLabels.end(); it++) {  
			m->mothurOut("Your file does not include the label " + *it); 
			if (processedLabels.count(lastLabel) != 1) {
				m->mothurOut(". I will use " + lastLabel + ".");  m->mothurOutEndLine();
				needToRun = true;
			}else {
				m->mothurOut(". Please refer to " + lastLabel + ".");  m->mothurOutEndLine();
			}
		}
		
		//run last label if you need to
		if (needToRun == true)  {
			if (sabund != NULL) {	delete sabund;	}
			sabund = (input->getSAbundVector(lastLabel));
			
			m->mothurOut(sabund->getLabel());  m->mothurOutEndLine();
			sabund->print(out);
			delete sabund;
			
			if (m->control_pressed) {  outputTypes.clear(); out.close(); remove(filename.c_str());  delete input; return 0; }
			
		}
		
		out.close();
		delete input;
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Name: "); m->mothurOutEndLine();
		m->mothurOut(filename); m->mothurOutEndLine();	outputNames.push_back(filename); outputTypes["sabund"].push_back(filename);
		m->mothurOutEndLine();
		
		//set sabund file as new current sabundfile
		string current = "";
		itTypes = outputTypes.find("sabund");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setSabundFile(current); }
		}
		
		return 0;		
	}

	catch(exception& e) {
		m->errorOut(e, "GetSAbundCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************


