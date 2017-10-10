/*
 *  getgroupcommand.cpp
 *  Mothur
 *
 *  Created by Thomas Ryabin on 2/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "getgroupcommand.h"
#include "inputdata.h"

//**********************************************************************************************************************
vector<string> GetgroupCommand::setParameters(){	
	try {
		CommandParameter pshared("shared", "InputTypes", "", "current", "none", "none", "none","",false,true, true); parameters.push_back(pshared);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "GetgroupCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string GetgroupCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The get.group command parameter is shared and it's required if you have no valid current file.\n";
		helpString += "You may not use any parameters with the get.group command.\n";
		helpString += "The get.group command should be in the following format: \n";
		helpString += "get.group()\n";
		helpString += "Example get.group().\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "GetgroupCommand", "getHelpString");
		exit(1);
	}
}

//**********************************************************************************************************************
GetgroupCommand::GetgroupCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
	}
	catch(exception& e) {
		m->errorOut(e, "GetgroupCommand", "GetgroupCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
GetgroupCommand::GetgroupCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		
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
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("shared");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["shared"] = inputDir + it->second;		}
				}
			}
			
			//get shared file
			sharedfile = validParameter.validFile(parameters, "shared", true);
			if (sharedfile == "not open") { sharedfile = ""; abort = true; }	
			else if (sharedfile == "not found") { 
				//if there is a current shared file, use it
				sharedfile = m->getSharedFile(); 
				if (sharedfile != "") { m->mothurOut("Using " + sharedfile + " as input file for the shared parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current sharedfile and the shared parameter is required."); m->mothurOutEndLine(); abort = true; }
			}else { m->setSharedFile(sharedfile); }
			
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = m->hasPath(sharedfile);		}
		}
	}
	catch(exception& e) {
		m->errorOut(e, "GetgroupCommand", "GetgroupCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int GetgroupCommand::execute(){
	try {
	
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
					
		InputData input(sharedfile, "sharedfile");
		SharedRAbundVectors* lookup = input.getSharedRAbundVectors();
        vector<string> namesOfGroups = lookup->getNamesGroups();
        delete lookup;
        
		for (int i = 0; i < namesOfGroups.size(); i++) {
			m->mothurOut(namesOfGroups[i]); m->mothurOutEndLine();
		}
    
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		m->mothurOutEndLine();
		
		return 0;	
	}

	catch(exception& e) {
		m->errorOut(e, "GetgroupCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************


