/*
 *  heatmapcommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 3/25/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "heatmapcommand.h"

//**********************************************************************************************************************
vector<string> HeatMapCommand::setParameters(){	
	try {
		CommandParameter plist("list", "InputTypes", "", "", "LRSS", "LRSS", "none","svg",false,false,true); parameters.push_back(plist);
		CommandParameter prabund("rabund", "InputTypes", "", "", "LRSS", "LRSS", "none","svg",false,false); parameters.push_back(prabund);
		CommandParameter psabund("sabund", "InputTypes", "", "", "LRSS", "LRSS", "none","svg",false,false); parameters.push_back(psabund);
		CommandParameter pshared("shared", "InputTypes", "", "", "LRSS", "LRSS", "none","svg",false,false,true); parameters.push_back(pshared);	
		CommandParameter prelabund("relabund", "InputTypes", "", "", "LRSS", "LRSS", "none","svg",false,false); parameters.push_back(prelabund);
		CommandParameter pgroups("groups", "String", "", "", "", "", "","",false,false); parameters.push_back(pgroups);
		CommandParameter pscale("scale", "Multiple", "log10-log2-linear", "log10", "", "", "","",false,false); parameters.push_back(pscale);
		CommandParameter psorted("sorted", "Multiple", "none-shared-topotu-topgroup", "shared", "", "", "","",false,false); parameters.push_back(psorted);
		CommandParameter pnumotu("numotu", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pnumotu);
		CommandParameter pfontsize("fontsize", "Number", "", "24", "", "", "","",false,false); parameters.push_back(pfontsize);
		CommandParameter plabel("label", "String", "", "", "", "", "","",false,false); parameters.push_back(plabel);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "HeatMapCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string HeatMapCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The heatmap.bin command parameters are shared, relabund, list, rabund, sabund, groups, sorted, scale, numotu, fontsize and label.  shared, relabund, list, rabund or sabund is required unless you have a valid current file.\n";
		helpString += "The groups parameter allows you to specify which of the groups in your groupfile you would like included in your heatmap.\n";
		helpString += "The sorted parameter allows you to order the otus displayed, default=shared, meaning display the shared otus first. Other options for sorted are none, meaning the exact representation of your otus, \n";
		helpString += "topotu, meaning the otus with the greatest abundance when totaled across groups, topgroup, meaning the top otus for each group. \n";
		helpString += "The scale parameter allows you to choose the range of color your bin information will be displayed with.\n";
		helpString += "The numotu parameter allows you to display only the top N otus, by default all the otus are displayed. You could choose to look at the top 10, by setting numotu=10. The default for sorted is topotu when numotu is used.\n";
		helpString += "The group names are separated by dashes. The label parameter allows you to select what distance levels you would like a heatmap created for, and are also separated by dashes.\n";
		helpString += "The fontsize parameter allows you to adjust the font size of the picture created, default=24.\n";
		helpString += "The heatmap.bin command should be in the following format: heatmap.bin(groups=yourGroups, sorted=yourSorted, label=yourLabels).\n";
		helpString += "Example heatmap.bin(groups=A-B-C, sorted=none, scale=log10).\n";
		helpString += "The default value for groups is all the groups in your groupfile, and all labels in your inputfile will be used.\n";
		helpString += "The default value for scale is log10; your other options are log2 and linear.\n";
		helpString += "The heatmap.bin command outputs a .svg file for each label you specify.\n";
		helpString += "Note: No spaces between parameter labels (i.e. groups), '=' and parameters (i.e.yourGroups).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "HeatMapCommand", "getHelpString");
		exit(1);
	}
}

//**********************************************************************************************************************
string HeatMapCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "svg") {  pattern = "[filename],svg"; } 
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "HeatMapCommand", "getOutputPattern");
        exit(1);
    }
}

//**********************************************************************************************************************
HeatMapCommand::HeatMapCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["svg"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "HeatMapCommand", "HeatMapCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

HeatMapCommand::HeatMapCommand(string option) {
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
			outputTypes["svg"] = tempOutNames;
		
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
				
				it = parameters.find("relabund");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["relabund"] = inputDir + it->second;		}
				}
			}
			
			//check for required parameters
			listfile = validParameter.validFile(parameters, "list", true);
			if (listfile == "not open") { abort = true; }
			else if (listfile == "not found") { listfile = ""; }
			else {  format = "list"; inputfile = listfile; m->setListFile(listfile); }
			
			sabundfile = validParameter.validFile(parameters, "sabund", true);
			if (sabundfile == "not open") {  abort = true; }	
			else if (sabundfile == "not found") { sabundfile = ""; }
			else {  format = "sabund"; inputfile = sabundfile; m->setSabundFile(sabundfile); }
			
			rabundfile = validParameter.validFile(parameters, "rabund", true);
			if (rabundfile == "not open") {  abort = true; }	
			else if (rabundfile == "not found") { rabundfile = ""; }
			else {  format = "rabund"; inputfile = rabundfile; m->setRabundFile(rabundfile); }
			
			sharedfile = validParameter.validFile(parameters, "shared", true);
			if (sharedfile == "not open") { abort = true; }	
			else if (sharedfile == "not found") { sharedfile = ""; }
			else {  format = "sharedfile"; inputfile = sharedfile; m->setSharedFile(sharedfile); }
			
			relabundfile = validParameter.validFile(parameters, "relabund", true);
			if (relabundfile == "not open") {  abort = true; }	
			else if (relabundfile == "not found") { relabundfile = ""; }
			else {  format = "relabund"; inputfile = relabundfile; m->setRelAbundFile(relabundfile); }
			
			
			if ((sharedfile == "") && (listfile == "") && (rabundfile == "") && (sabundfile == "") && (relabundfile == "")) { 
				//is there are current file available for any of these?
				//give priority to shared, then list, then rabund, then sabund
				//if there is a current shared file, use it
				sharedfile = m->getSharedFile(); 
				if (sharedfile != "") { inputfile = sharedfile; format = "sharedfile"; m->mothurOut("Using " + sharedfile + " as input file for the shared parameter."); m->mothurOutEndLine(); }
				else { 
					listfile = m->getListFile(); 
					if (listfile != "") { inputfile = listfile; format = "list"; m->mothurOut("Using " + listfile + " as input file for the list parameter."); m->mothurOutEndLine(); }
					else { 
						rabundfile = m->getRabundFile(); 
						if (rabundfile != "") { inputfile = rabundfile; format = "rabund"; m->mothurOut("Using " + rabundfile + " as input file for the rabund parameter."); m->mothurOutEndLine(); }
						else { 
							sabundfile = m->getSabundFile(); 
							if (sabundfile != "") { inputfile = sabundfile; format = "sabund"; m->mothurOut("Using " + sabundfile + " as input file for the sabund parameter."); m->mothurOutEndLine(); }
							else { 
								relabundfile = m->getRelAbundFile(); 
								if (relabundfile != "") { inputfile = relabundfile; format = "relabund"; m->mothurOut("Using " + relabundfile + " as input file for the relabund parameter."); m->mothurOutEndLine(); }
								else { 
									m->mothurOut("No valid current files. You must provide a list, sabund, rabund, relabund or shared file."); m->mothurOutEndLine(); 
									abort = true;
								}
							}
						}
					}
				}
			}
			
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = m->hasPath(inputfile);		}
			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  m->splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
			
			
			groups = validParameter.validFile(parameters, "groups", false);			
			if (groups == "not found") { groups = ""; }
			else { 
				m->splitAtDash(groups, Groups);
                if (Groups.size() != 0) { if (Groups[0] != "all") { Groups.clear(); } }
			}
			
			string temp = validParameter.validFile(parameters, "numotu", false);		if (temp == "not found") { temp = "0"; }
			m->mothurConvert(temp, numOTU);
			
			temp = validParameter.validFile(parameters, "fontsize", false);				if (temp == "not found") { temp = "24"; }
			m->mothurConvert(temp, fontSize);
			
			sorted = validParameter.validFile(parameters, "sorted", false);				
			if (sorted == "not found") { 
				//if numOTU is used change default
				if (numOTU != 0) { sorted = "topotu"; }
				else { sorted = "shared"; }
			}
		 
			scale = validParameter.validFile(parameters, "scale", false);				if (scale == "not found") { scale = "log10"; }
			
			if ((sorted != "none") && (sorted != "shared") && (sorted != "topotu") && (sorted != "topgroup")) { m->mothurOut(sorted + " is not a valid sorting option. Sorted options are: none, shared, topotu, topgroup"); m->mothurOutEndLine(); abort=true;  }
		}

	}
	catch(exception& e) {
		m->errorOut(e, "HeatMapCommand", "HeatMapCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int HeatMapCommand::execute(){
	try {
	
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
		string lastLabel;
		input = new InputData(inputfile, format, Groups);
		
        vector<string> currentLabels;
		if (format == "sharedfile") {
			//you have groups
			lookup = input->getSharedRAbundVectors();
			lastLabel = lookup->getLabel();
            currentLabels = lookup->getOTUNames();
            Groups = lookup->getNamesGroups();
		}else if ((format == "list") || (format == "rabund") || (format == "sabund")) {
			//you are using just a list file and have only one group
			rabund = input->getRAbundVector();
			lastLabel = rabund->getLabel();
		}else if (format == "relabund") {
			//you have groups
			lookupFloat = input->getSharedRAbundFloatVectors();
			lastLabel = lookupFloat->getLabel();
            currentLabels = lookup->getOTUNames();
            Groups = lookupFloat->getNamesGroups();
		}
		
        heatmap = new HeatMap(sorted, scale, numOTU, fontSize, outputDir, inputfile, currentLabels);
        
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;

		if (format == "sharedfile") {	
		
			//as long as you are not at the end of the file or done wih the lines you want
			while((lookup != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
				if (m->getControl_pressed()) {
                    delete lookup;
					for (int i = 0; i < outputNames.size(); i++) {	if (outputNames[i] != "control") {  m->mothurRemove(outputNames[i]);  } } outputTypes.clear();
					 
					delete input; delete heatmap; return 0;
				}
		
				if(allLines == 1 || labels.count(lookup->getLabel()) == 1){
	
					m->mothurOut(lookup->getLabel()); m->mothurOutEndLine();
                    vector<SharedRAbundVector*> data = lookup->getSharedRAbundVectors();
					string outputFileName = heatmap->getPic(data, lookup->getNamesGroups());
                    for (int i = 0; i < data.size(); i++) {  delete data[i];  }
					outputNames.push_back(outputFileName); outputTypes["svg"].push_back(outputFileName);
					
					processedLabels.insert(lookup->getLabel());
					userLabels.erase(lookup->getLabel());
				}
				
				if ((m->anyLabelsToProcess(lookup->getLabel(), userLabels, "") ) && (processedLabels.count(lastLabel) != 1)) {
					string saveLabel = lookup->getLabel();
                    delete lookup;
			
					lookup = input->getSharedRAbundVectors(lastLabel);
					m->mothurOut(lookup->getLabel()); m->mothurOutEndLine();
					
                    vector<SharedRAbundVector*> data = lookup->getSharedRAbundVectors();
					string outputFileName = heatmap->getPic(data, lookup->getNamesGroups());
                    for (int i = 0; i < data.size(); i++) {  delete data[i];  }
					outputNames.push_back(outputFileName); outputTypes["svg"].push_back(outputFileName);
					
					processedLabels.insert(lookup->getLabel());
					userLabels.erase(lookup->getLabel());
					
					//restore real lastlabel to save below
					lookup->setLabels(saveLabel);
				}
				
				lastLabel = lookup->getLabel();
				//prevent memory leak
				delete lookup;
							
				//get next line to process
				lookup = input->getSharedRAbundVectors();
			}
			
			
			if (m->getControl_pressed()) {
				for (int i = 0; i < outputNames.size(); i++) {	if (outputNames[i] != "control") {  m->mothurRemove(outputNames[i]);  } } outputTypes.clear();
				 
				delete input; delete heatmap; return 0;
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
			if (needToRun )  {
				delete lookup;
				lookup = input->getSharedRAbundVectors(lastLabel);
				
				m->mothurOut(lookup->getLabel()); m->mothurOutEndLine();
                vector<SharedRAbundVector*> data = lookup->getSharedRAbundVectors();
                string outputFileName = heatmap->getPic(data, lookup->getNamesGroups());
                for (int i = 0; i < data.size(); i++) {  delete data[i];  }

				outputNames.push_back(outputFileName); outputTypes["svg"].push_back(outputFileName);
				delete lookup;
			}
		
			//reset groups parameter
			  
			
		}else if ((format == "list") || (format == "rabund") || (format == "sabund")) {
	
			while((rabund != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
				if (m->getControl_pressed()) {   
					for (int i = 0; i < outputNames.size(); i++) {	if (outputNames[i] != "control") {  m->mothurRemove(outputNames[i]);  } } outputTypes.clear();
					delete rabund;  delete input; delete heatmap; return 0;	
				}

				if(allLines == 1 || labels.count(rabund->getLabel()) == 1){			
	
					m->mothurOut(rabund->getLabel()); m->mothurOutEndLine();
					string outputFileName = heatmap->getPic(rabund);
					outputNames.push_back(outputFileName); outputTypes["svg"].push_back(outputFileName);
					
					processedLabels.insert(rabund->getLabel());
					userLabels.erase(rabund->getLabel());
				}
				
				if ((m->anyLabelsToProcess(rabund->getLabel(), userLabels, "") ) && (processedLabels.count(lastLabel) != 1)) {
					string saveLabel = rabund->getLabel();
					
					delete rabund;
					rabund = input->getRAbundVector(lastLabel);
					m->mothurOut(rabund->getLabel()); m->mothurOutEndLine();
					
					string outputFileName = heatmap->getPic(rabund);
					outputNames.push_back(outputFileName); outputTypes["svg"].push_back(outputFileName);
					
					processedLabels.insert(rabund->getLabel());
					userLabels.erase(rabund->getLabel());
					
					//restore real lastlabel to save below
					rabund->setLabel(saveLabel);
				}		
				
								
								
				lastLabel = rabund->getLabel();			
				delete rabund;
				rabund = input->getRAbundVector();
			}
			
			if (m->getControl_pressed()) {
				for (int i = 0; i < outputNames.size(); i++) {	if (outputNames[i] != "control") {  m->mothurRemove(outputNames[i]);  } } outputTypes.clear();
				delete input; delete heatmap; return 0;
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
			if (needToRun )  {
		
				if (rabund != NULL) {	delete rabund;	}
				rabund = input->getRAbundVector(lastLabel);
				m->mothurOut(rabund->getLabel()); m->mothurOutEndLine();
					
				string outputFileName = heatmap->getPic(rabund);
				outputNames.push_back(outputFileName); outputTypes["svg"].push_back(outputFileName);
				delete rabund; 
			}
		
		}else {
		
			//as long as you are not at the end of the file or done wih the lines you want
			while((lookupFloat != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
				if (m->getControl_pressed()) {
                    delete lookupFloat;
					for (int i = 0; i < outputNames.size(); i++) {	if (outputNames[i] != "control") {  m->mothurRemove(outputNames[i]);  } } outputTypes.clear();
					 
					delete input; delete heatmap; return 0;
				}
		
				if(allLines == 1 || labels.count(lookupFloat->getLabel()) == 1){
	
					m->mothurOut(lookupFloat->getLabel()); m->mothurOutEndLine();
                    vector<SharedRAbundFloatVector*> data = lookupFloat->getSharedRAbundFloatVectors();
                    string outputFileName = heatmap->getPic(data, lookupFloat->getNamesGroups());
                    for (int i = 0; i < data.size(); i++) {  delete data[i];  }
					outputNames.push_back(outputFileName); outputTypes["svg"].push_back(outputFileName);
					
					processedLabels.insert(lookupFloat->getLabel());
					userLabels.erase(lookupFloat->getLabel());
				}
				
				if ((m->anyLabelsToProcess(lookupFloat->getLabel(), userLabels, "") ) && (processedLabels.count(lastLabel) != 1)) {
					string saveLabel = lookupFloat->getLabel();
				
					delete lookupFloat;
					lookupFloat = input->getSharedRAbundFloatVectors(lastLabel);
					m->mothurOut(lookupFloat->getLabel()); m->mothurOutEndLine();
					
                    vector<SharedRAbundFloatVector*> data = lookupFloat->getSharedRAbundFloatVectors();
                    string outputFileName = heatmap->getPic(data, lookupFloat->getNamesGroups());
                    for (int i = 0; i < data.size(); i++) {  delete data[i];  }
					outputNames.push_back(outputFileName); outputTypes["svg"].push_back(outputFileName);
					
					processedLabels.insert(lookupFloat->getLabel());
					userLabels.erase(lookupFloat->getLabel());
					
					//restore real lastlabel to save below
					lookupFloat->setLabels(saveLabel);
				}
				
				lastLabel = lookupFloat->getLabel();
				//prevent memory leak
				delete lookupFloat;
							
				//get next line to process
				lookupFloat = input->getSharedRAbundFloatVectors();
			}
			
			
			if (m->getControl_pressed()) {
				for (int i = 0; i < outputNames.size(); i++) {	if (outputNames[i] != "control") {  m->mothurRemove(outputNames[i]);  } } outputTypes.clear();
				 
				delete input; delete heatmap; return 0;
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
			if (needToRun )  {
                delete lookupFloat;
                lookupFloat = input->getSharedRAbundFloatVectors(lastLabel);
                m->mothurOut(lookupFloat->getLabel()); m->mothurOutEndLine();
                
                vector<SharedRAbundFloatVector*> data = lookupFloat->getSharedRAbundFloatVectors();
                string outputFileName = heatmap->getPic(data, lookupFloat->getNamesGroups());
                for (int i = 0; i < data.size(); i++) {  delete data[i];  }

				outputNames.push_back(outputFileName); outputTypes["svg"].push_back(outputFileName);
				delete lookupFloat;
			}
		
        }
		
		delete input; 
		delete heatmap;
		
		if (m->getControl_pressed()) {
			for (int i = 0; i < outputNames.size(); i++) {	if (outputNames[i] != "control") {  m->mothurRemove(outputNames[i]);  } } outputTypes.clear(); return 0;
		}
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();

		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "HeatMapCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************


