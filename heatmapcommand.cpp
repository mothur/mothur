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

HeatMapCommand::HeatMapCommand(string option) {
	try {
		globaldata = GlobalData::getInstance();
		abort = false;
		allLines = 1;
		labels.clear();
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string AlignArray[] =  {"groups","label","sorted","scale","fontsize","numotu","outputdir","inputdir"};
			vector<string> myArray (AlignArray, AlignArray+(sizeof(AlignArray)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			
			//check to make sure all parameters are valid for command
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += m->hasPath(globaldata->inputFileName); //if user entered a file with a path then preserve it	
			}
			
			//make sure the user has already run the read.otu command
			if ((globaldata->getListFile() == "") && (globaldata->getRabundFile() == "") && (globaldata->getSabundFile() == "") && (globaldata->getSharedFile() == "") && (globaldata->getRelAbundFile() == "")) {
				 m->mothurOut("You must read a list, rabund, sabund, or a list and a group, shared, or relabund file before you can use the heatmap.bin command."); m->mothurOutEndLine(); abort = true; 
			}

			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  m->splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
			
			//if the user has not specified any labels use the ones from read.otu
			if (label == "") {  
				allLines = globaldata->allLines; 
				labels = globaldata->labels; 
			}
			
			groups = validParameter.validFile(parameters, "groups", false);			
			if (groups == "not found") { groups = ""; }
			else { 
				m->splitAtDash(groups, Groups);
				globaldata->Groups = Groups;
			}
			
			string temp = validParameter.validFile(parameters, "numotu", false);		if (temp == "not found") { temp = "0"; }
			convert(temp, numOTU);
			
			temp = validParameter.validFile(parameters, "fontsize", false);				if (temp == "not found") { temp = "24"; }
			convert(temp, fontSize);
			
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

void HeatMapCommand::help(){
	try {
		m->mothurOut("The heatmap.bin command can only be executed after a successful read.otu command.\n");
		m->mothurOut("The heatmap.bin command parameters are groups, sorted, scale, numotu, fontsize and label.  No parameters are required.\n");
		m->mothurOut("The groups parameter allows you to specify which of the groups in your groupfile you would like included in your heatmap.\n");
		m->mothurOut("The sorted parameter allows you to order the otus displayed, default=shared, meaning display the shared otus first. Other options for sorted are none, meaning the exact representation of your otus, \n");
		m->mothurOut("topotu, meaning the otus with the greatest abundance when totaled across groups, topgroup, meaning the top otus for each group. \n");
		m->mothurOut("The scale parameter allows you to choose the range of color your bin information will be displayed with.\n");
		m->mothurOut("The numotu parameter allows you to display only the top N otus, by default all the otus are displayed. You could choose to look at the top 10, by setting numotu=10. The default for sorted is topotu when numotu is used.\n");
		m->mothurOut("The group names are separated by dashes. The label parameter allows you to select what distance levels you would like a heatmap created for, and are also separated by dashes.\n");
		m->mothurOut("The fontsize parameter allows you to adjust the font size of the picture created, default=24.\n");
		m->mothurOut("The heatmap.bin command should be in the following format: heatmap.bin(groups=yourGroups, sorted=yourSorted, label=yourLabels).\n");
		m->mothurOut("Example heatmap.bin(groups=A-B-C, sorted=none, scale=log10).\n");
		m->mothurOut("The default value for groups is all the groups in your groupfile, and all labels in your inputfile will be used.\n");
		m->mothurOut("The default value for scale is log10; your other options are log2 and linear.\n");
		m->mothurOut("The heatmap.bin command outputs a .svg file for each label you specify.\n");
		m->mothurOut("Note: No spaces between parameter labels (i.e. groups), '=' and parameters (i.e.yourGroups).\n\n");

	}
	catch(exception& e) {
		m->errorOut(e, "HeatMapCommand", "help");
		exit(1);
	}
}

//**********************************************************************************************************************

HeatMapCommand::~HeatMapCommand(){
}

//**********************************************************************************************************************

int HeatMapCommand::execute(){
	try {
	
		if (abort == true) { return 0; }
		
		heatmap = new HeatMap(sorted, scale, numOTU, fontSize, outputDir);
		format = globaldata->getFormat();

		string lastLabel;
		
		read = new ReadOTUFile(globaldata->inputFileName);	
		read->read(&*globaldata); 
		input = globaldata->ginput;
		
		if (format == "sharedfile") {
			//you have groups
			lookup = input->getSharedRAbundVectors();
			lastLabel = lookup[0]->getLabel();
	
		}else if ((format == "list") || (format == "rabund") || (format == "sabund")) {
			//you are using just a list file and have only one group
			rabund = globaldata->rabund;
			lastLabel = rabund->getLabel();
		}else if (format == "relabund") {
			//you have groups
			lookupFloat = input->getSharedRAbundFloatVectors();
			lastLabel = lookupFloat[0]->getLabel();
		}
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;

		if (format == "sharedfile") {	
		
			//as long as you are not at the end of the file or done wih the lines you want
			while((lookup[0] != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
				if (m->control_pressed) {
					for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }
					for (int i = 0; i < outputNames.size(); i++) {	if (outputNames[i] != "control") {  remove(outputNames[i].c_str());  } }
					globaldata->Groups.clear(); 
					delete read; delete heatmap; return 0;
				}
		
				if(allLines == 1 || labels.count(lookup[0]->getLabel()) == 1){			
	
					m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
					outputNames.push_back(heatmap->getPic(lookup));
					
					processedLabels.insert(lookup[0]->getLabel());
					userLabels.erase(lookup[0]->getLabel());
				}
				
				if ((m->anyLabelsToProcess(lookup[0]->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
					string saveLabel = lookup[0]->getLabel();
				
					for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }  
					lookup = input->getSharedRAbundVectors(lastLabel);
					m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
					
					outputNames.push_back(heatmap->getPic(lookup));
					
					processedLabels.insert(lookup[0]->getLabel());
					userLabels.erase(lookup[0]->getLabel());
					
					//restore real lastlabel to save below
					lookup[0]->setLabel(saveLabel);
				}
				
				lastLabel = lookup[0]->getLabel();
				//prevent memory leak
				for (int i = 0; i < lookup.size(); i++) {  delete lookup[i]; lookup[i] = NULL; }
							
				//get next line to process
				lookup = input->getSharedRAbundVectors();				
			}
			
			
			if (m->control_pressed) {
				for (int i = 0; i < outputNames.size(); i++) {	if (outputNames[i] != "control") {  remove(outputNames[i].c_str());  } }
				globaldata->Groups.clear(); 
				delete read; delete heatmap; return 0;
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
			if (needToRun == true)  {
				for (int i = 0; i < lookup.size(); i++) { if (lookup[i] != NULL) { delete lookup[i]; } }  
				lookup = input->getSharedRAbundVectors(lastLabel);
				
				m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
				outputNames.push_back(heatmap->getPic(lookup));
				for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }
			}
		
			//reset groups parameter
			globaldata->Groups.clear();  
			
		}else if ((format == "list") || (format == "rabund") || (format == "sabund")) {
	
			while((rabund != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
				if (m->control_pressed) {   
					for (int i = 0; i < outputNames.size(); i++) {	if (outputNames[i] != "control") {  remove(outputNames[i].c_str());  } }
					delete rabund;  delete read; delete heatmap; return 0;	
				}

				if(allLines == 1 || labels.count(rabund->getLabel()) == 1){			
	
					m->mothurOut(rabund->getLabel()); m->mothurOutEndLine();
					outputNames.push_back(heatmap->getPic(rabund));
					
					processedLabels.insert(rabund->getLabel());
					userLabels.erase(rabund->getLabel());
				}
				
				if ((m->anyLabelsToProcess(rabund->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
					string saveLabel = rabund->getLabel();
					
					delete rabund;
					rabund = input->getRAbundVector(lastLabel);
					m->mothurOut(rabund->getLabel()); m->mothurOutEndLine();
					
					outputNames.push_back(heatmap->getPic(rabund));
					
					processedLabels.insert(rabund->getLabel());
					userLabels.erase(rabund->getLabel());
					
					//restore real lastlabel to save below
					rabund->setLabel(saveLabel);
				}		
				
								
								
				lastLabel = rabund->getLabel();			
				delete rabund;
				rabund = input->getRAbundVector();
			}
			
			if (m->control_pressed) {
				for (int i = 0; i < outputNames.size(); i++) {	if (outputNames[i] != "control") {  remove(outputNames[i].c_str());  } }
				delete read; delete heatmap; return 0;
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
			if (needToRun == true)  {
		
				if (rabund != NULL) {	delete rabund;	}
				rabund = input->getRAbundVector(lastLabel);
				m->mothurOut(rabund->getLabel()); m->mothurOutEndLine();
					
				outputNames.push_back(heatmap->getPic(rabund));
				delete rabund; globaldata->rabund = NULL;
			}
		
		}else {
		
			//as long as you are not at the end of the file or done wih the lines you want
			while((lookupFloat[0] != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
				if (m->control_pressed) {
					for (int i = 0; i < lookupFloat.size(); i++) {  delete lookupFloat[i];  }
					for (int i = 0; i < outputNames.size(); i++) {	if (outputNames[i] != "control") {  remove(outputNames[i].c_str());  } }
					globaldata->Groups.clear(); 
					delete read; delete heatmap; return 0;
				}
		
				if(allLines == 1 || labels.count(lookupFloat[0]->getLabel()) == 1){			
	
					m->mothurOut(lookupFloat[0]->getLabel()); m->mothurOutEndLine();
					outputNames.push_back(heatmap->getPic(lookupFloat));
					
					processedLabels.insert(lookupFloat[0]->getLabel());
					userLabels.erase(lookupFloat[0]->getLabel());
				}
				
				if ((m->anyLabelsToProcess(lookupFloat[0]->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
					string saveLabel = lookupFloat[0]->getLabel();
				
					for (int i = 0; i < lookupFloat.size(); i++) {  delete lookupFloat[i];  }  
					lookupFloat = input->getSharedRAbundFloatVectors(lastLabel);
					m->mothurOut(lookupFloat[0]->getLabel()); m->mothurOutEndLine();
					
					outputNames.push_back(heatmap->getPic(lookupFloat));
					
					processedLabels.insert(lookupFloat[0]->getLabel());
					userLabels.erase(lookupFloat[0]->getLabel());
					
					//restore real lastlabel to save below
					lookupFloat[0]->setLabel(saveLabel);
				}
				
				lastLabel = lookupFloat[0]->getLabel();
				//prevent memory leak
				for (int i = 0; i < lookupFloat.size(); i++) {  delete lookupFloat[i]; lookupFloat[i] = NULL; }
							
				//get next line to process
				lookupFloat = input->getSharedRAbundFloatVectors();				
			}
			
			
			if (m->control_pressed) {
				for (int i = 0; i < outputNames.size(); i++) {	if (outputNames[i] != "control") {  remove(outputNames[i].c_str());  } }
				globaldata->Groups.clear(); 
				delete read; delete heatmap; return 0;
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
			if (needToRun == true)  {
				for (int i = 0; i < lookupFloat.size(); i++) { if (lookupFloat[i] != NULL) { delete lookupFloat[i]; } }  
				lookupFloat = input->getSharedRAbundFloatVectors(lastLabel);
				
				m->mothurOut(lookupFloat[0]->getLabel()); m->mothurOutEndLine();
				outputNames.push_back(heatmap->getPic(lookupFloat));
				for (int i = 0; i < lookupFloat.size(); i++) {  delete lookupFloat[i];  }
			}
		
			//reset groups parameter
			globaldata->Groups.clear();  

		}
		
		globaldata->rabund = NULL;
		delete input; globaldata->ginput = NULL;
		
		if (m->control_pressed) {
			for (int i = 0; i < outputNames.size(); i++) {	if (outputNames[i] != "control") {  remove(outputNames[i].c_str());  } }
			delete read; delete heatmap; return 0;
		}
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
	
		delete read;
		delete heatmap;

		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "HeatMapCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************


