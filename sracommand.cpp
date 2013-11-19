//
//  sracommand.cpp
//  Mothur
//
//  Created by SarahsWork on 10/28/13.
//  Copyright (c) 2013 Schloss Lab. All rights reserved.
//

#include "sracommand.h"

//**********************************************************************************************************************
vector<string> SRACommand::setParameters(){
	try {
        CommandParameter psff("sff", "InputTypes", "", "", "none", "none", "none","sra",false,false); parameters.push_back(psff);
		CommandParameter pfastqfile("fastqfile", "InputTypes", "", "", "none", "none", "none","sra",false,false); parameters.push_back(pfastqfile);
        //choose only one multiple options
        CommandParameter pplatform("platform", "Multiple", "454-???-???", "454", "", "", "","",false,false); parameters.push_back(pplatform);
         //every command must have inputdir and outputdir.  This allows mothur users to redirect input and output files.
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "SRACommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string SRACommand::getHelpString(){
	try {
		string helpString = "";
		helpString += "The sra command creates a sequence read archive from sff or fastq files.\n";
		helpString += "The sra command parameters are: sff, fastqfiles, oligos, platform....\n";
		helpString += "The sffiles parameter is used to provide a file containing a \n";
		helpString += "The new command should be in the following format: \n";
		helpString += "new(...)\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "SRACommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string SRACommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "sra") {  pattern = "[filename],sra"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "SRACommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
SRACommand::SRACommand(){
	try {
		abort = true; calledHelp = true;
		setParameters();
        vector<string> tempOutNames;
		outputTypes["sra"] = tempOutNames; 
	}
	catch(exception& e) {
		m->errorOut(e, "SRACommand", "SRACommand");
		exit(1);
	}
}
//**********************************************************************************************************************
SRACommand::SRACommand(string option)  {
	try {
		abort = false; calledHelp = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			//valid paramters for this command
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			map<string,string>::iterator it;
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) {
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			
			//if the user changes the input directory command factory will send this info to us in the output parameter
			string inputDir = validParameter.validFile(parameters, "inputdir", false);
			if (inputDir == "not found"){	inputDir = "";		}
			else {
            
                string path;
				it = parameters.find("sfffiles");
				//user has given a template file
				if(it != parameters.end()){
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["sfffiles"] = inputDir + it->second;		}
				}
				
				it = parameters.find("fastqfiles");
				//user has given a template file
				if(it != parameters.end()){
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["fastqfiles"] = inputDir + it->second;		}
				}
            }
            
			//check for parameters
            fastqfiles = validParameter.validFile(parameters, "fastqfiles", true);
			if (fastqfiles == "not open") { fastqfiles = "";  abort = true; }
			else if (fastqfiles == "not found") { fastqfiles = ""; }
			
			sfffiles = validParameter.validFile(parameters, "sfffiles", true);
			if (sfffiles == "not open") {  sfffiles = "";  abort = true; }
			else if (sfffiles == "not found") { sfffiles = ""; }
			
			if ((fastqfiles == "") && (sfffiles == "")) {
                m->mothurOut("No valid current files. You must provide a sfffiles or fastqfiles file before you can use the sra command."); m->mothurOutEndLine(); abort = true;
            }
			            
            //use only one Mutliple type
			platform = validParameter.validFile(parameters, "platform", false);
			if (platform == "not found") { platform = "454"; }
			
			if ((platform == "454") || (platform == "????") || (platform == "????") || (platform == "????")) { }
			else { m->mothurOut("Not a valid platform option.  Valid platform options are 454, ...."); m->mothurOutEndLine(); abort = true; }
            
            			
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "SRACommand", "SRACommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int SRACommand::execute(){
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
        

		
        //output files created by command
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
        return 0;
		
    }
	catch(exception& e) {
		m->errorOut(e, "SRACommand", "SRACommand");
		exit(1);
	}
}
//**********************************************************************************************************************


