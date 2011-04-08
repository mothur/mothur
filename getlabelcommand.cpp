/*
 *  GetlabelCommand.cpp
 *  Mothur
 *
 *  Created by Thomas Ryabin on 1/30/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "getlabelcommand.h"


//**********************************************************************************************************************
vector<string> GetlabelCommand::setParameters(){	
	try {
		CommandParameter plist("list", "InputTypes", "", "", "LRSS", "LRSS", "none",false,false); parameters.push_back(plist);
		CommandParameter prabund("rabund", "InputTypes", "", "", "LRSS", "LRSS", "none",false,false); parameters.push_back(prabund);
		CommandParameter psabund("sabund", "InputTypes", "", "", "LRSS", "LRSS", "none",false,false); parameters.push_back(psabund);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "GetlabelCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
GetlabelCommand::GetlabelCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
	}
	catch(exception& e) {
		m->errorOut(e, "GetlabelCommand", "CollectCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
string GetlabelCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The get.label command parameters are list, sabund and rabund file. \n";
		helpString += "The get.label command should be in the following format: \n";
		helpString += "get.label()\n";
		helpString += "Example get.label().\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "GetlabelCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************

GetlabelCommand::GetlabelCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		
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
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				
				it = parameters.find("rabund");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["rabund"] = inputDir + it->second;		}
				}
				
				it = parameters.find("sabund");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["sabund"] = inputDir + it->second;		}
				}
				
				it = parameters.find("list");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["list"] = inputDir + it->second;		}
				}
			}
			
			//check for required parameters
			listfile = validParameter.validFile(parameters, "list", true);
			if (listfile == "not open") { listfile = ""; abort = true; }
			else if (listfile == "not found") { listfile = ""; }
			else {  format = "list"; inputfile = listfile; }
			
			sabundfile = validParameter.validFile(parameters, "sabund", true);
			if (sabundfile == "not open") { sabundfile = ""; abort = true; }	
			else if (sabundfile == "not found") { sabundfile = ""; }
			else {  format = "sabund"; inputfile = sabundfile; }
			
			rabundfile = validParameter.validFile(parameters, "rabund", true);
			if (rabundfile == "not open") { rabundfile = ""; abort = true; }	
			else if (rabundfile == "not found") { rabundfile = ""; }
			else {  format = "rabund"; inputfile = rabundfile; }
			
			if ((listfile == "") && (rabundfile == "") && (sabundfile == "")) { 
				//is there are current file available for any of these?
				//give priority to list, then rabund, then sabund
				//if there is a current shared file, use it
				
				listfile = m->getListFile(); 
				if (listfile != "") { inputfile = listfile; format = "list"; m->mothurOut("Using " + listfile + " as input file for the list parameter."); m->mothurOutEndLine(); }
				else { 
					rabundfile = m->getRabundFile(); 
					if (rabundfile != "") { inputfile = rabundfile; format = "rabund"; m->mothurOut("Using " + rabundfile + " as input file for the rabund parameter."); m->mothurOutEndLine(); }
					else { 
						sabundfile = m->getSabundFile(); 
						if (sabundfile != "") { inputfile = sabundfile; format = "sabund"; m->mothurOut("Using " + sabundfile + " as input file for the sabund parameter."); m->mothurOutEndLine(); }
						else { 
							m->mothurOut("No valid current files. You must provide a list, sabund or rabund file."); m->mothurOutEndLine(); 
							abort = true;
						}
					}
				}
			}
		}

	}
	catch(exception& e) {
		m->errorOut(e, "GetlabelCommand", "GetlabelCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int GetlabelCommand::execute(){
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		InputData* input = new InputData(inputfile, format);
		OrderVector* order = input->getOrderVector();
		string label = order->getLabel();
		
		while (order != NULL) {
			
			if (m->control_pressed) { delete input;  delete order; return 0; }
			
			label = order->getLabel();	
			
			m->mothurOut(label); m->mothurOutEndLine();
			
			delete order;		
			order = input->getOrderVector();
		}
		
		delete input; 
		
		return 0;	
	}

	catch(exception& e) {
		m->errorOut(e, "GetlabelCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************


