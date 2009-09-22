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

HeatMapCommand::HeatMapCommand(string option){
	try {
		globaldata = GlobalData::getInstance();
		abort = false;
		allLines = 1;
		labels.clear();
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string AlignArray[] =  {"groups","label","sorted","scale"};
			vector<string> myArray (AlignArray, AlignArray+(sizeof(AlignArray)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			
			//check to make sure all parameters are valid for command
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//make sure the user has already run the read.otu command
			if ((globaldata->getListFile() == "") && (globaldata->getRabundFile() == "") && (globaldata->getSabundFile() == "") && (globaldata->getSharedFile() == "")) {
				 mothurOut("You must read a list, rabund, sabund, or a list and a group, or a shared before you can use the heatmap.bin command."); mothurOutEndLine(); abort = true; 
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
			
			if (abort == false) {
				heatmap = new HeatMap(sorted, scale);
				format = globaldata->getFormat();
			}
		}

	}
	catch(exception& e) {
		errorOut(e, "HeatMapCommand", "HeatMapCommand");
		exit(1);
	}
}

//**********************************************************************************************************************

void HeatMapCommand::help(){
	try {
		mothurOut("The heatmap.bin command can only be executed after a successful read.otu command.\n");
		mothurOut("The heatmap.bin command parameters are groups, sorted, scale label.  No parameters are required.\n");
		mothurOut("The groups parameter allows you to specify which of the groups in your groupfile you would like included in your heatmap.\n");
		mothurOut("The sorted parameter allows you to choose to see the file with the shared otus at the top or the otus in the order they appear in your input file. \n");
		mothurOut("The scale parameter allows you to choose the range of color your bin information will be displayed with.\n");
		mothurOut("The group names are separated by dashes. The label parameter allows you to select what distance levels you would like a heatmap created for, and are also separated by dashes.\n");
		mothurOut("The heatmap.bin command should be in the following format: heatmap.bin(groups=yourGroups, sorted=yourSorted, label=yourLabels).\n");
		mothurOut("Example heatmap.bin(groups=A-B-C, sorted=F, scale=log10).\n");
		mothurOut("The default value for groups is all the groups in your groupfile, and all labels in your inputfile will be used.\n");
		mothurOut("The default value for sorted is T meaning you want the shared otus on top, you may change it to F meaning the exact representation of your input file.\n");
		mothurOut("The default value for scale is log10; your other options are log2 and linear.\n");
		mothurOut("The heatmap.bin command outputs a .svg file for each label you specify.\n");
		mothurOut("Note: No spaces between parameter labels (i.e. groups), '=' and parameters (i.e.yourGroups).\n\n");

	}
	catch(exception& e) {
		errorOut(e, "HeatMapCommand", "help");
		exit(1);
	}
}

//**********************************************************************************************************************

HeatMapCommand::~HeatMapCommand(){
	if (abort == false) {
		delete read;
		delete heatmap;
	}
}

//**********************************************************************************************************************

int HeatMapCommand::execute(){
	try {
	
		if (abort == true) { return 0; }

		string lastLabel;
	
		if (format == "sharedfile") {
			//you have groups
			read = new ReadOTUFile(globaldata->inputFileName);	
			read->read(&*globaldata); 
			
			input = globaldata->ginput;
			lookup = input->getSharedRAbundVectors();
			lastLabel = lookup[0]->getLabel();
	
		}else if ((format == "list") || (format == "rabund") || (format == "sabund")) {
			//you are using just a list file and have only one group
			read = new ReadOTUFile(globaldata->inputFileName);	
			read->read(&*globaldata); 
		
			rabund = globaldata->rabund;
			lastLabel = rabund->getLabel();
			input = globaldata->ginput;	

		}
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;

		if ((format != "list") && (format != "rabund") && (format != "sabund")) {	
		
			//as long as you are not at the end of the file or done wih the lines you want
			while((lookup[0] != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
		
				if(allLines == 1 || labels.count(lookup[0]->getLabel()) == 1){			
	
					mothurOut(lookup[0]->getLabel()); mothurOutEndLine();
					heatmap->getPic(lookup);
					
					processedLabels.insert(lookup[0]->getLabel());
					userLabels.erase(lookup[0]->getLabel());
				}
				
				if ((anyLabelsToProcess(lookup[0]->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
					for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }  
					lookup = input->getSharedRAbundVectors(lastLabel);
					mothurOut(lookup[0]->getLabel()); mothurOutEndLine();
					
					heatmap->getPic(lookup);
					
					processedLabels.insert(lookup[0]->getLabel());
					userLabels.erase(lookup[0]->getLabel());
				}
				
				lastLabel = lookup[0]->getLabel();
				//prevent memory leak
				for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }
							
				//get next line to process
				lookup = input->getSharedRAbundVectors();				
			}
			
			//output error messages about any remaining user labels
			set<string>::iterator it;
			bool needToRun = false;
			for (it = userLabels.begin(); it != userLabels.end(); it++) {  
				mothurOut("Your file does not include the label " + *it); 
				if (processedLabels.count(lastLabel) != 1) {
					mothurOut(". I will use " + lastLabel + "."); mothurOutEndLine();
					needToRun = true;
				}else {
					mothurOut(". Please refer to " + lastLabel + "."); mothurOutEndLine();
				}
			}
		
			//run last label if you need to
			if (needToRun == true)  {
				for (int i = 0; i < lookup.size(); i++) { if (lookup[i] != NULL) { delete lookup[i]; } }  
				lookup = input->getSharedRAbundVectors(lastLabel);
				
				mothurOut(lookup[0]->getLabel()); mothurOutEndLine();
				heatmap->getPic(lookup);
				for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }
			}
		
			
			
			//reset groups parameter
			globaldata->Groups.clear();  
			
		}else{
	
			while((rabund != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {

				if(allLines == 1 || labels.count(rabund->getLabel()) == 1){			
	
					mothurOut(rabund->getLabel()); mothurOutEndLine();
					heatmap->getPic(rabund);
					
					processedLabels.insert(rabund->getLabel());
					userLabels.erase(rabund->getLabel());
				}
				
				if ((anyLabelsToProcess(rabund->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {

					delete rabund;
					rabund = input->getRAbundVector(lastLabel);
					mothurOut(rabund->getLabel()); mothurOutEndLine();
					
					heatmap->getPic(rabund);
					
					processedLabels.insert(rabund->getLabel());
					userLabels.erase(rabund->getLabel());
				}		
				
								
								
				lastLabel = rabund->getLabel();			
				delete rabund;
				rabund = input->getRAbundVector();
			}
			
			//output error messages about any remaining user labels
			set<string>::iterator it;
			bool needToRun = false;
			for (it = userLabels.begin(); it != userLabels.end(); it++) {  
				mothurOut("Your file does not include the label " + *it); 
				if (processedLabels.count(lastLabel) != 1) {
					mothurOut(". I will use " + lastLabel + "."); mothurOutEndLine();
					needToRun = true;
				}else {
					mothurOut(". Please refer to " + lastLabel + "."); mothurOutEndLine();
				}
			}
		
			//run last label if you need to
			if (needToRun == true)  {
		
				if (rabund != NULL) {	delete rabund;	}
				rabund = input->getRAbundVector(lastLabel);
				mothurOut(rabund->getLabel()); mothurOutEndLine();
					
				heatmap->getPic(rabund);
				delete rabund; globaldata->rabund = NULL;
			}
		
		}
		
		globaldata->rabund = NULL;
		delete input; globaldata->ginput = NULL;
		return 0;
	}
	catch(exception& e) {
		errorOut(e, "HeatMapCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************


