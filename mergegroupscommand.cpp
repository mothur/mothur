/*
 *  mergegroupscommand.cpp
 *  mothur
 *
 *  Created by westcott on 1/24/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "mergegroupscommand.h"
#include "sharedutilities.h"

//**********************************************************************************************************************
vector<string> MergeGroupsCommand::setParameters(){	
	try {
		CommandParameter pshared("shared", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(pshared);
		CommandParameter pdesign("design", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(pdesign);
		CommandParameter plabel("label", "String", "", "", "", "", "",false,false); parameters.push_back(plabel);
		CommandParameter pgroups("groups", "String", "", "", "", "", "",false,false); parameters.push_back(pgroups);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "MergeGroupsCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string MergeGroupsCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The merge.groups command reads a shared file and a design file and merges the groups in the shared file that are in the same grouping in the design file.\n";
		helpString += "The merge.groups command outputs a .shared file. \n";
		helpString += "The merge.groups command parameters are shared, groups, label and design.  The design and shared parameter are required.\n";
		helpString += "The design parameter allows you to assign your groups to sets. It is required. \n";
		helpString += "The design file looks like the group file.  It is a 2 column tab delimited file, where the first column is the group name and the second column is the set the group belongs to.\n";
		helpString += "The groups parameter allows you to specify which of the groups in your shared you would like included. The group names are separated by dashes.\n";
		helpString += "The label parameter allows you to select what distance levels you would like, and are also separated by dashes.\n";
		helpString += "The merge.groups command should be in the following format: merge.groups(design=yourDesignFile, shared=yourSharedFile).\n";
		helpString += "Example merge.groups(design=temp.design, groups=A-B-C, shared=temp.shared).\n";
		helpString += "The default value for groups is all the groups in your sharedfile, and all labels in your inputfile will be used.\n";
		helpString += "Note: No spaces between parameter labels (i.e. groups), '=' and parameters (i.e.yourGroups).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "MergeGroupsCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
MergeGroupsCommand::MergeGroupsCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["shared"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "MergeGroupsCommand", "MetaStatsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

MergeGroupsCommand::MergeGroupsCommand(string option) {
	try {
		abort = false; calledHelp = false;   
		allLines = 1;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			
			//check to make sure all parameters are valid for command
			map<string,string>::iterator it;
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["shared"] = tempOutNames;
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";	}
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("design");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["design"] = inputDir + it->second;		}
				}
				
				it = parameters.find("shared");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["shared"] = inputDir + it->second;		}
				}
			}
			
			//check for required parameters
			designfile = validParameter.validFile(parameters, "design", true);
			if (designfile == "not open") { abort = true; }
			else if (designfile == "not found") {  				
				//if there is a current shared file, use it
				designfile = m->getDesignFile(); 
				if (designfile != "") { m->mothurOut("Using " + designfile + " as input file for the design parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current designfile and the design parameter is required."); m->mothurOutEndLine(); abort = true; }
			}	
			
			//make sure the user has already run the read.otu command
			sharedfile = validParameter.validFile(parameters, "shared", true);
			if (sharedfile == "not open") { abort = true; }
			else if (sharedfile == "not found") {  				
				//if there is a current shared file, use it
				sharedfile = m->getSharedFile(); 
				if (sharedfile != "") { m->mothurOut("Using " + sharedfile + " as input file for the shared parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current sharedfile and the shared parameter is required."); m->mothurOutEndLine(); abort = true; }
			}	
			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  m->splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
			
			groups = validParameter.validFile(parameters, "groups", false);			
			if (groups == "not found") { groups = "all";  }
			m->splitAtDash(groups, Groups);
			m->Groups = Groups;
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "MergeGroupsCommand", "MergeGroupsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int MergeGroupsCommand::execute(){
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		if (outputDir == "") {  outputDir += m->hasPath(sharedfile);  }
		string outputFileName = outputDir + m->getRootName(m->getSimpleName(sharedfile)) + "merge" +  m->getExtension(sharedfile);
		outputTypes["shared"].push_back(outputFileName); outputNames.push_back(outputFileName);
		
		ofstream out;
		m->openOutputFile(outputFileName, out);
		
		designMap = new GroupMap(designfile);
		designMap->readDesignMap();
		
		InputData input(sharedfile, "sharedfile");
		lookup = input.getSharedRAbundVectors();
		string lastLabel = lookup[0]->getLabel();
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;
		
		//as long as you are not at the end of the file or done wih the lines you want
		while((lookup[0] != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			
			if (m->control_pressed) {  out.close(); for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } m->Groups.clear();  delete designMap;  for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); } return 0; }
			
			if(allLines == 1 || labels.count(lookup[0]->getLabel()) == 1){			
				
				m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
				process(lookup, out);
				
				processedLabels.insert(lookup[0]->getLabel());
				userLabels.erase(lookup[0]->getLabel());
			}
			
			if ((m->anyLabelsToProcess(lookup[0]->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
				string saveLabel = lookup[0]->getLabel();
				
				for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }  
				lookup = input.getSharedRAbundVectors(lastLabel);
				m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
				
				process(lookup, out);
				
				processedLabels.insert(lookup[0]->getLabel());
				userLabels.erase(lookup[0]->getLabel());
				
				//restore real lastlabel to save below
				lookup[0]->setLabel(saveLabel);
			}
			
			lastLabel = lookup[0]->getLabel();
			//prevent memory leak
			for (int i = 0; i < lookup.size(); i++) {  delete lookup[i]; lookup[i] = NULL; }
			
			if (m->control_pressed) {  out.close(); m->Groups.clear();   delete designMap;  for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); } return 0; }
			
			//get next line to process
			lookup = input.getSharedRAbundVectors();				
		}
		
		if (m->control_pressed) { out.close(); m->Groups.clear();  delete designMap;  for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); }  return 0; }
		
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
			for (int i = 0; i < lookup.size(); i++) { if (lookup[i] != NULL) { delete lookup[i]; } }  
			lookup = input.getSharedRAbundVectors(lastLabel);
			
			m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
			
			process(lookup, out);
			
			for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }
		}
		
		out.close();
		//reset groups parameter
		m->Groups.clear();  
		delete designMap;
		
		if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); } return 0;}
		
		
		//set shared file as new current sharedfile
		string current = "";
		itTypes = outputTypes.find("shared");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setSharedFile(current); }
		}
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "MergeGroupsCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************

int MergeGroupsCommand::process(vector<SharedRAbundVector*>& thisLookUp, ofstream& out){
	try {
		
		map<string, SharedRAbundVector> merged;
		map<string, SharedRAbundVector>::iterator it;
		
		for (int i = 0; i < thisLookUp.size(); i++) {
			
			if (m->control_pressed) { return 0; }
			
			//what grouping does this group belong to
			string grouping = designMap->getGroup(thisLookUp[i]->getGroup());
			if (grouping == "not found") { m->mothurOut("[ERROR]: " + thisLookUp[i]->getGroup() + " is not in your design file. Ignoring!"); m->mothurOutEndLine(); grouping = "NOTFOUND"; }
			
			else {
				//do we already have a member of this grouping?
				it = merged.find(grouping);
				
				if (it == merged.end()) { //nope, so create it
					merged[grouping] = *thisLookUp[i];
					merged[grouping].setGroup(grouping);
				}else { //yes, merge it
					
					for (int j = 0; j < thisLookUp[i]->getNumBins(); j++) {
						int abund = (it->second).getAbundance(j);
						abund += thisLookUp[i]->getAbundance(j);
						
						(it->second).set(j, abund, grouping);
					}
				}
			}
		}
		
		//print new file
		for (it = merged.begin(); it != merged.end(); it++) {
			out << (it->second).getLabel() << '\t' << it->first << '\t';
			(it->second).print(out);
		}
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "MergeGroupsCommand", "process");
		exit(1);
	}
}
//**********************************************************************************************************************



