/*
 *  getrelabundcommand.cpp
 *  Mothur
 *
 *  Created by westcott on 6/21/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "getrelabundcommand.h"

//**********************************************************************************************************************
vector<string> GetRelAbundCommand::setParameters(){	
	try {
		CommandParameter pshared("shared", "InputTypes", "", "", "none", "none", "none","relabund",false,true, true); parameters.push_back(pshared);
		CommandParameter pgroups("groups", "String", "", "", "", "", "","",false,false); parameters.push_back(pgroups);
		CommandParameter pscale("scale", "Multiple", "totalgroup-totalotu-averagegroup-averageotu", "totalgroup", "", "", "","",false,false); parameters.push_back(pscale);
		CommandParameter plabel("label", "String", "", "", "", "", "","",false,false); parameters.push_back(plabel);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "GetRelAbundCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string GetRelAbundCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The get.relabund command parameters are shared, groups, scale and label.  shared is required, unless you have a valid current file.\n";
		helpString += "The groups parameter allows you to specify which of the groups in your groupfile you would like included. The group names are separated by dashes.\n";
		helpString += "The label parameter allows you to select what distance levels you would like, and are also separated by dashes.\n";
		helpString += "The scale parameter allows you to select what scale you would like to use. Choices are totalgroup, totalotu, averagegroup, averageotu, default is totalgroup.\n";
		helpString += "The get.relabund command should be in the following format: get.relabund(groups=yourGroups, label=yourLabels).\n";
		helpString += "Example get.relabund(groups=A-B-C, scale=averagegroup).\n";
		helpString += "The default value for groups is all the groups in your groupfile, and all labels in your inputfile will be used.\n";
		helpString += "The get.relabund command outputs a .relabund file.\n";
		
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "GetRelAbundCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string GetRelAbundCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "relabund")      {   pattern = "[filename],relabund";    }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "GetRelAbundCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
GetRelAbundCommand::GetRelAbundCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["relabund"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "GetRelAbundCommand", "GetRelAbundCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
GetRelAbundCommand::GetRelAbundCommand(string option) {
	try {
		abort = false; calledHelp = false;   
		allLines = true;
				
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
				if (!validParameter.isValidParameter(it->first, myArray, it->second)) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["relabund"] = tempOutNames;
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.valid(parameters, "inputdir");		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("shared");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["shared"] = inputDir + it->second;		}
				}
			}
		
			//get shared file
			sharedfile = validParameter.validFile(parameters, "shared");
			if (sharedfile == "not open") { sharedfile = ""; abort = true; }	
			else if (sharedfile == "not found") { 
				//if there is a current shared file, use it
				sharedfile = current->getSharedFile(); 
				if (sharedfile != "") { m->mothurOut("Using " + sharedfile + " as input file for the shared parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current sharedfile and the shared parameter is required."); m->mothurOutEndLine(); abort = true; }
			}else { current->setSharedFile(sharedfile); }
			
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){	outputDir = util.hasPath(sharedfile);		}

			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			label = validParameter.valid(parameters, "label");			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  util.splitAtDash(label, labels);  allLines = false;  }
				else { allLines = true;  }
			}
			
			groups = validParameter.valid(parameters, "groups");			
			if (groups == "not found") { groups = ""; pickedGroups = false; }
			else { 
				pickedGroups = true;
				util.splitAtDash(groups, Groups);
                    if (Groups.size() != 0) { if (Groups[0]== "all") { Groups.clear(); } }
			}
			
			scale = validParameter.valid(parameters, "scale");				if (scale == "not found") { scale = "totalgroup"; }
			
			if ((scale != "totalgroup") && (scale != "totalotu") && (scale != "averagegroup") && (scale != "averageotu")) {
				m->mothurOut(scale + " is not a valid scaling option for the get.relabund command. Choices are totalgroup, totalotu, averagegroup, averageotu."); m->mothurOutEndLine(); abort = true; 
			}
		}

	}
	catch(exception& e) {
		m->errorOut(e, "GetRelAbundCommand", "GetRelAbundCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int GetRelAbundCommand::execute(){
	try {
	
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
        map<string, string> variables; 
		variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(sharedfile));
		string outputFileName = getOutputFileName("relabund", variables);
		ofstream out;
		util.openOutputFile(outputFileName, out);
		out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);
		
        InputData input(sharedfile, "sharedfile", Groups);
        set<string> processedLabels;
        set<string> userLabels = labels;
        string lastLabel = "";
        
        SharedRAbundVectors* lookup = util.getNextShared(input, allLines, userLabels, processedLabels, lastLabel);
        Groups = lookup->getNamesGroups();

        vector<string> binLabels = lookup->getOTUNames();
        out << "label\tGroup\tnumOtus";
        for (int i = 0; i < binLabels.size(); i++) { out  << '\t' << binLabels[i]; } out << endl;
        
        while (lookup != NULL) {
            
            if (m->getControl_pressed()) { delete lookup; break; }
            
            getRelAbundance(lookup, out); delete lookup;
            
            lookup = util.getNextShared(input, allLines, userLabels, processedLabels, lastLabel);
        }
        out.close();
			
		if (m->getControl_pressed()) { outputTypes.clear(); util.mothurRemove(outputFileName); return 0;    }
		
		m->mothurOut("\nOutput File Names: \n"); 
		m->mothurOut(outputFileName); m->mothurOutEndLine(); outputNames.push_back(outputFileName); outputTypes["relabund"].push_back(outputFileName);
		m->mothurOutEndLine();
		
		//set relabund file as new current relabundfile
		string currentName = "";
		itTypes = outputTypes.find("relabund");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setRelAbundFile(currentName); }
		}
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "GetRelAbundCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************

int GetRelAbundCommand::getRelAbundance(SharedRAbundVectors*& thisLookUp, ofstream& out){
	try {
        vector<string> groups = thisLookUp->getNamesGroups();
        vector<string> binLabels = thisLookUp->getOTUNames();
        
		 for (int i = 0; i < thisLookUp->size(); i++) {
			out << thisLookUp->getLabel() << '\t' << groups[i] << '\t' << thisLookUp->getNumBins();
			
			for (int j = 0; j < thisLookUp->getNumBins(); j++) {
			
				if (m->getControl_pressed()) { return 0; }
			
				int abund = thisLookUp->get(j, groups[i]);
				
				float relabund = 0.0;
				
				if (scale == "totalgroup") { 
					relabund = abund / (float) thisLookUp->getNumSeqs(groups[i]);
				}else if (scale == "totalotu") {
					//calc the total in this otu
					int totalOtu = thisLookUp->getOTUTotal(j);
					relabund = abund / (float) totalOtu;
				}else if (scale == "averagegroup") {
					relabund = abund / (float) (thisLookUp->getNumSeqs(groups[i]) / (float) thisLookUp->getNumBins());
				}else if (scale == "averageotu") {
					//calc the total in this otu
					int totalOtu = thisLookUp->getOTUTotal(j);
					float averageOtu = totalOtu / (float) thisLookUp->size();
					
					relabund = abund / (float) averageOtu;
				}else{ m->mothurOut(scale + " is not a valid scaling option."); m->mothurOutEndLine(); m->setControl_pressed(true); return 0; }
				
				out  << '\t' << relabund;
			}
			out << endl;
		 }
	
		 return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "GetRelAbundCommand", "getRelAbundance");
		exit(1);
	}
}
//**********************************************************************************************************************


