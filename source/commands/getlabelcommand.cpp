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
		CommandParameter plist("list", "InputTypes", "", "", "LRSS", "LRSS", "none","",false,false, true); parameters.push_back(plist);
		CommandParameter prabund("rabund", "InputTypes", "", "", "LRSS", "LRSS", "none","",false,false, true); parameters.push_back(prabund);
		CommandParameter psabund("sabund", "InputTypes", "", "", "LRSS", "LRSS", "none","",false,false, true); parameters.push_back(psabund);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
        
        abort = false; calledHelp = false;
		
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

GetlabelCommand::GetlabelCommand(string option) : Command()  {
	try {
    
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }
		
		else {
			OptionParser parser(option, setParameters());
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			listfile = validParameter.validFile(parameters, "list");
			if (listfile == "not open") { listfile = ""; abort = true; }
			else if (listfile == "not found") { listfile = ""; }
			else {  format = "list"; inputfile = listfile; current->setListFile(listfile); }
			
			sabundfile = validParameter.validFile(parameters, "sabund");
			if (sabundfile == "not open") { sabundfile = ""; abort = true; }	
			else if (sabundfile == "not found") { sabundfile = ""; }
			else {  format = "sabund"; inputfile = sabundfile; current->setSabundFile(sabundfile); }
			
			rabundfile = validParameter.validFile(parameters, "rabund");
			if (rabundfile == "not open") { rabundfile = ""; abort = true; }	
			else if (rabundfile == "not found") { rabundfile = ""; }
			else {  format = "rabund"; inputfile = rabundfile; current->setRabundFile(rabundfile); }
			
			if ((listfile == "") && (rabundfile == "") && (sabundfile == "")) { 
				//is there are current file available for any of these?
				//give priority to list, then rabund, then sabund
				//if there is a current shared file, use it
				
				listfile = current->getListFile(); 
				if (listfile != "") { inputfile = listfile; format = "list"; m->mothurOut("Using " + listfile + " as input file for the list parameter.\n");  }
				else { 
					rabundfile = current->getRabundFile(); 
					if (rabundfile != "") { inputfile = rabundfile; format = "rabund"; m->mothurOut("Using " + rabundfile + " as input file for the rabund parameter.\n");  }
					else { 
						sabundfile = current->getSabundFile(); 
						if (sabundfile != "") { inputfile = sabundfile; format = "sabund"; m->mothurOut("Using " + sabundfile + " as input file for the sabund parameter.\n");  }
						else { 
							m->mothurOut("No valid current files. You must provide a list, sabund or rabund file.\n");  
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
		
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
		InputData* input = new InputData(inputfile, format, nullVector);
		OrderVector* order = input->getOrderVector();
		string label = order->getLabel();
		
		while (order != nullptr) {
			
			if (m->getControl_pressed()) { delete input;  delete order; return 0; }
			
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


