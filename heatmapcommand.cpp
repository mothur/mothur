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
			string AlignArray[] =  {"groups","label","sorted","scale","outputdir","inputdir"};
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
				outputDir += hasPath(globaldata->inputFileName); //if user entered a file with a path then preserve it	
			}
			
			//make sure the user has already run the read.otu command
			if ((globaldata->getListFile() == "") && (globaldata->getRabundFile() == "") && (globaldata->getSabundFile() == "") && (globaldata->getSharedFile() == "")) {
				 m->mothurOut("You must read a list, rabund, sabund, or a list and a group, or a shared before you can use the heatmap.bin command."); m->mothurOutEndLine(); abort = true; 
			}

			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  splitAtDash(label, labels);  allLines = 0;  }
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
				splitAtDash(groups, Groups);
				globaldata->Groups = Groups;
			}
			
			sorted = validParameter.validFile(parameters, "sorted", false);			if (sorted == "not found") { sorted = "T"; }
		 
			scale = validParameter.validFile(parameters, "scale", false);				if (scale == "not found") { scale = "log10"; }
			
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
		m->mothurOut("The heatmap.bin command parameters are groups, sorted, scale label.  No parameters are required.\n");
		m->mothurOut("The groups parameter allows you to specify which of the groups in your groupfile you would like included in your heatmap.\n");
		m->mothurOut("The sorted parameter allows you to choose to see the file with the shared otus at the top or the otus in the order they appear in your input file. \n");
		m->mothurOut("The scale parameter allows you to choose the range of color your bin information will be displayed with.\n");
		m->mothurOut("The group names are separated by dashes. The label parameter allows you to select what distance levels you would like a heatmap created for, and are also separated by dashes.\n");
		m->mothurOut("The heatmap.bin command should be in the following format: heatmap.bin(groups=yourGroups, sorted=yourSorted, label=yourLabels).\n");
		m->mothurOut("Example heatmap.bin(groups=A-B-C, sorted=F, scale=log10).\n");
		m->mothurOut("The default value for groups is all the groups in your groupfile, and all labels in your inputfile will be used.\n");
		m->mothurOut("The default value for sorted is T meaning you want the shared otus on top, you may change it to F meaning the exact representation of your input file.\n");
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
		
		heatmap = new HeatMap(sorted, scale, outputDir);
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
		}
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;

		if ((format != "list") && (format != "rabund") && (format != "sabund")) {	
		
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
				
				if ((anyLabelsToProcess(lookup[0]->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
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
			
		}else{
	
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
				
				if ((anyLabelsToProcess(rabund->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
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


