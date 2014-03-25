//
//  getmimarkspackagecommand.cpp
//  Mothur
//
//  Created by Sarah Westcott on 3/25/14.
//  Copyright (c) 2014 Schloss Lab. All rights reserved.
//

#include "getmimarkspackagecommand.h"

//**********************************************************************************************************************
vector<string> GetMIMarksPackageCommand::setParameters(){
	try {
        //files that have dependancies
        CommandParameter pgroup("group", "InputTypes", "", "", "groupOligos", "none", "none","",false,false); parameters.push_back(pgroup);
        CommandParameter poligos("oligos", "InputTypes", "", "", "groupOligos", "none", "none","",false,false); parameters.push_back(poligos);
  		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "GetMIMarksPackageCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string GetMIMarksPackageCommand::getHelpString(){
	try {
		string helpString = "";
		helpString += "The get.mimarkspackage command creates a mimarks package form with your groups. The required fields are flagged with * characters.\n";
        helpString += "Further documentation on the different packages and required formats can be found here, http://www.mothur.org/wiki/MIMarks_Data_Packages.\n";
		helpString += "The get.mimarkspackage command parameters are: oligos, group and package. oligos or group is required.\n";
		helpString += "The oligos parameter is used to provide your oligos file so mothur can extract your group names.\n";
        helpString += "The group parameter is used to provide your group file so mothur can extract your group names.\n";
        helpString += "The package parameter is used to select the mimarks package you would like to use. Default=???\n";
		helpString += "The get.mimarkspackage command should be in the following format: get.mimarkspackage(oligos=yourOligosFile, package=yourPackage)\n";
		helpString += "get.mimarkspackage(oligos=GQY1XT001.oligos, package=???)\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "GetMIMarksPackageCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string GetMIMarksPackageCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "tsv") {  pattern = "[filename],tsv"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "GetMIMarksPackageCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
GetMIMarksPackageCommand::GetMIMarksPackageCommand(){
	try {
		abort = true; calledHelp = true;
		setParameters();
        vector<string> tempOutNames;
		outputTypes["tsv"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "GetMIMarksPackageCommand", "GetMIMarksPackageCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
GetMIMarksPackageCommand::GetMIMarksPackageCommand(string option)  {
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
			
            vector<string> tempOutNames;
			outputTypes["tsv"] = tempOutNames;
            
			//if the user changes the input directory command factory will send this info to us in the output parameter
			string inputDir = validParameter.validFile(parameters, "inputdir", false);
			if (inputDir == "not found"){	inputDir = "";		}
			else {
                
				string path;
				it = parameters.find("oligos");
				//user has given a template file
				if(it != parameters.end()){
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["oligos"] = inputDir + it->second;		}
				}
				
				it = parameters.find("group");
				//user has given a template file
				if(it != parameters.end()){
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["group"] = inputDir + it->second;		}
				}
				
            }
            
			groupfile = validParameter.validFile(parameters, "group", true);
			if (groupfile == "not open") {  groupfile = "";  abort = true; }
			else if (groupfile == "not found") { groupfile = ""; }
            else {  m->setGroupFile(groupfile); }
            
            oligosfile = validParameter.validFile(parameters, "oligos", true);
			if (oligosfile == "not found")      {	oligosfile = "";	}
			else if(oligosfile == "not open")	{	abort = true;		}
			else {	m->setOligosFile(oligosfile); }

            if ((groupfile != "") && (oligosfile != "")) {
                m->mothurOut("[ERROR]: You may not use a group file and an oligos file, only one."); m->mothurOutEndLine(); abort = true;
            }

            if ((groupfile == "") && (oligosfile == "")) {
                oligosfile = m->getOligosFile();
                if (oligosfile != "") {  m->mothurOut("Using " + oligosfile + " as input file for the oligos parameter."); m->mothurOutEndLine(); }
                else {
                    groupfile = m->getGroupFile();
                    if (groupfile != "") {  m->mothurOut("Using " + groupfile + " as input file for the group parameter."); m->mothurOutEndLine(); }
                    else {
                        m->mothurOut("[ERROR]: You must provide groupfile or oligos file for the get.mimarkspackage command."); m->mothurOutEndLine(); abort = true;
                    }
                }
            }
            
            package = validParameter.validFile(parameters, "package", false);         if (package == "not found") { package = "package"; }
            //if (!checkCasesPackage(package)) { abort = true; } //error message in function
            
            //turn _ to spaces mothur's work around
            for (int i = 0; i < package.length(); i++) { if (package[i] == '_') { package[i] = ' '; }  }


		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "GetMIMarksPackageCommand", "GetMIMarksPackageCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int GetMIMarksPackageCommand::execute(){
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
		m->errorOut(e, "GetMIMarksPackageCommand", "GetMIMarksPackageCommand");
		exit(1);
	}
}
//**********************************************************************************************************************


