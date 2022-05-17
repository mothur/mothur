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
        
        abort = false; calledHelp = false; allLines = true;
        
        vector<string> tempOutNames;
        outputTypes["svg"] = tempOutNames;
		
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

HeatMapCommand::HeatMapCommand(string option) : Command() {
	try {

		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }
		
		else {
			OptionParser parser(option, setParameters());
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			listfile = validParameter.validFile(parameters, "list");
			if (listfile == "not open") { abort = true; }
			else if (listfile == "not found") { listfile = ""; }
			else {  format = "list"; inputfile = listfile; current->setListFile(listfile); }
			
			sabundfile = validParameter.validFile(parameters, "sabund");
			if (sabundfile == "not open") {  abort = true; }	
			else if (sabundfile == "not found") { sabundfile = ""; }
			else {  format = "sabund"; inputfile = sabundfile; current->setSabundFile(sabundfile); }
			
			rabundfile = validParameter.validFile(parameters, "rabund");
			if (rabundfile == "not open") {  abort = true; }	
			else if (rabundfile == "not found") { rabundfile = ""; }
			else {  format = "rabund"; inputfile = rabundfile; current->setRabundFile(rabundfile); }
			
			sharedfile = validParameter.validFile(parameters, "shared");
			if (sharedfile == "not open") { abort = true; }	
			else if (sharedfile == "not found") { sharedfile = ""; }
			else {  format = "sharedfile"; inputfile = sharedfile; current->setSharedFile(sharedfile); }
			
			relabundfile = validParameter.validFile(parameters, "relabund");
			if (relabundfile == "not open") {  abort = true; }	
			else if (relabundfile == "not found") { relabundfile = ""; }
			else {  format = "relabund"; inputfile = relabundfile; current->setRelAbundFile(relabundfile); }
			
			
			if ((sharedfile == "") && (listfile == "") && (rabundfile == "") && (sabundfile == "") && (relabundfile == "")) { 
				//is there are current file available for any of these?
				//give priority to shared, then list, then rabund, then sabund
				//if there is a current shared file, use it
				sharedfile = current->getSharedFile(); 
				if (sharedfile != "") { inputfile = sharedfile; format = "sharedfile"; m->mothurOut("Using " + sharedfile + " as input file for the shared parameter.\n");  }
				else { 
					listfile = current->getListFile(); 
					if (listfile != "") { inputfile = listfile; format = "list"; m->mothurOut("Using " + listfile + " as input file for the list parameter.\n");  }
					else { 
						rabundfile = current->getRabundFile(); 
						if (rabundfile != "") { inputfile = rabundfile; format = "rabund"; m->mothurOut("Using " + rabundfile + " as input file for the rabund parameter.\n");  }
						else { 
							sabundfile = current->getSabundFile(); 
							if (sabundfile != "") { inputfile = sabundfile; format = "sabund"; m->mothurOut("Using " + sabundfile + " as input file for the sabund parameter.\n");  }
							else { 
								relabundfile = current->getRelAbundFile(); 
								if (relabundfile != "") { inputfile = relabundfile; format = "relabund"; m->mothurOut("Using " + relabundfile + " as input file for the relabund parameter.\n");  }
								else { 
									m->mothurOut("No valid current files. You must provide a list, sabund, rabund, relabund or shared file.\n");  
									abort = true;
								}
							}
						}
					}
				}
			}
			
			 
					if (outputdir == ""){    outputdir = util.hasPath(inputfile);		}
			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			label = validParameter.valid(parameters, "label");			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  util.splitAtDash(label, labels);  allLines = false;  }
				else { allLines = true;  }
			}
			
			groups = validParameter.valid(parameters, "groups");			
			if (groups == "not found") { groups = ""; }
			else { 
				util.splitAtDash(groups, Groups);
                if (Groups.size() != 0) { if (Groups[0]== "all") { Groups.clear(); } }
			}
			
			string temp = validParameter.valid(parameters, "numotu");		if (temp == "not found") { temp = "0"; }
			util.mothurConvert(temp, numOTU);
			
			temp = validParameter.valid(parameters, "fontsize");				if (temp == "not found") { temp = "24"; }
			util.mothurConvert(temp, fontSize);
			
			sorted = validParameter.valid(parameters, "sorted");
			if (sorted == "not found") { 
				//if numOTU is used change default
				if (numOTU != 0) { sorted = "topotu"; }
				else { sorted = "shared"; }
			}
		 
			scale = validParameter.valid(parameters, "scale");				if (scale == "not found") { scale = "log10"; }
			
			if ((sorted != "none") && (sorted != "shared") && (sorted != "topotu") && (sorted != "topgroup")) { m->mothurOut(sorted + " is not a valid sorting option. Sorted options are: none, shared, topotu, topgroup\n");  abort=true;  }
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
        if (abort) { if (calledHelp) { return 0; }  return 2;    }
        
        InputData input(inputfile, format, Groups);
        set<string> processedLabels;
        set<string> userLabels = labels;
        string lastLabel = "";
        
        HeatMap heatmap(sorted, scale, numOTU, fontSize, outputdir, inputfile);
		
		if (format == "sharedfile") {
            SharedRAbundVectors* lookup = util.getNextShared(input, allLines, userLabels, processedLabels, lastLabel);
            Groups = lookup->getNamesGroups();
        
            while (lookup != nullptr) {
                
                if (m->getControl_pressed()) { delete lookup; break; }
                
                string outputFileName = heatmap.getPic(lookup); delete lookup;
                outputNames.push_back(outputFileName); outputTypes["svg"].push_back(outputFileName);
                
                lookup = util.getNextShared(input, allLines, userLabels, processedLabels, lastLabel);
            }
            
		}else if ((format == "list") || (format == "rabund") || (format == "sabund")) {
			RAbundVector* rabund = util.getNextRAbund(input, allLines, userLabels, processedLabels, lastLabel);
                   
            while (rabund != nullptr) {
                       
                if (m->getControl_pressed()) { delete rabund; break; }
                       
                string outputFileName = heatmap.getPic(rabund); delete rabund;
                outputNames.push_back(outputFileName); outputTypes["svg"].push_back(outputFileName);
                      
                rabund = util.getNextRAbund(input, allLines, userLabels, processedLabels, lastLabel);
            }
            
		}else if (format == "relabund") {
            SharedRAbundFloatVectors* lookup = util.getNextRelabund(input, allLines, userLabels, processedLabels, lastLabel);
            Groups = lookup->getNamesGroups();
            
            while (lookup != nullptr) {
                
                if (m->getControl_pressed()) { delete lookup; break; }
                
                string outputFileName = heatmap.getPic(lookup); delete lookup;
                outputNames.push_back(outputFileName); outputTypes["svg"].push_back(outputFileName);
                
                lookup = util.getNextRelabund(input, allLines, userLabels, processedLabels, lastLabel);
            }
		}

		if (m->getControl_pressed()) {
			for (int i = 0; i < outputNames.size(); i++) {	if (outputNames[i] != "control") {  util.mothurRemove(outputNames[i]);  } } outputTypes.clear(); return 0;
		}
		
		m->mothurOut("\nOutput File Names: \n"); 
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i] +"\n"); 	} m->mothurOutEndLine();

		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "HeatMapCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************


