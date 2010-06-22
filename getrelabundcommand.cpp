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

GetRelAbundCommand::GetRelAbundCommand(string option) {
	try {
		globaldata = GlobalData::getInstance();
		abort = false;
		allLines = 1;
		labels.clear();
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string AlignArray[] =  {"groups","label","scale","outputdir","inputdir"};
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
			if ((globaldata->getSharedFile() == "")) {
				 m->mothurOut("You must read a list and a group, or a shared file before you can use the get.relabund command."); m->mothurOutEndLine(); abort = true; 
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
			
			scale = validParameter.validFile(parameters, "scale", false);				if (scale == "not found") { scale = "totalgroup"; }
			
			if ((scale != "totalgroup") && (scale != "totalgroup") && (scale != "totalgroup") && (scale != "totalgroup")) {
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

void GetRelAbundCommand::help(){
	try {
		m->mothurOut("The get.relabund command can only be executed after a successful read.otu command of a list and group or shared file.\n");
		m->mothurOut("The get.relabund command parameters are groups, scale and label.  No parameters are required.\n");
		m->mothurOut("The groups parameter allows you to specify which of the groups in your groupfile you would like included. The group names are separated by dashes.\n");
		m->mothurOut("The label parameter allows you to select what distance levels you would like, and are also separated by dashes.\n");
		m->mothurOut("The scale parameter allows you to select what scale you would like to use. Choices are totalgroup, totalotu, averagegroup, averageotu, default is totalgroup.\n");
		m->mothurOut("The get.relabund command should be in the following format: get.relabund(groups=yourGroups, label=yourLabels).\n");
		m->mothurOut("Example get.relabund(groups=A-B-C, scale=log10).\n");
		m->mothurOut("The default value for groups is all the groups in your groupfile, and all labels in your inputfile will be used.\n");
		m->mothurOut("The get.relabund command outputs a .relabund file.\n");
		m->mothurOut("Note: No spaces between parameter labels (i.e. groups), '=' and parameters (i.e.yourGroups).\n\n");

	}
	catch(exception& e) {
		m->errorOut(e, "GetRelAbundCommand", "help");
		exit(1);
	}
}

//**********************************************************************************************************************

GetRelAbundCommand::~GetRelAbundCommand(){
}

//**********************************************************************************************************************

int GetRelAbundCommand::execute(){
	try {
	
		if (abort == true) { return 0; }
		
		string outputFileName = outputDir + getRootName(getSimpleName(globaldata->inputFileName)) + "relabund";
		ofstream out;
		openOutputFile(outputFileName, out);
		
		read = new ReadOTUFile(globaldata->inputFileName);	
		read->read(&*globaldata); 
		input = globaldata->ginput;
		lookup = input->getSharedRAbundVectors();
		string lastLabel = lookup[0]->getLabel();
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;

		//as long as you are not at the end of the file or done wih the lines you want
		while((lookup[0] != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			
			if (m->control_pressed) {  for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } globaldata->Groups.clear(); delete read;  out.close(); remove(outputFileName.c_str()); return 0; }
	
			if(allLines == 1 || labels.count(lookup[0]->getLabel()) == 1){			

				m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
				getRelAbundance(lookup, out);
				
				processedLabels.insert(lookup[0]->getLabel());
				userLabels.erase(lookup[0]->getLabel());
			}
			
			if ((anyLabelsToProcess(lookup[0]->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
				string saveLabel = lookup[0]->getLabel();
			
				for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }  
				lookup = input->getSharedRAbundVectors(lastLabel);
				m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
				
				getRelAbundance(lookup, out);
				
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
		
		if (m->control_pressed) { globaldata->Groups.clear(); delete read;  out.close(); remove(outputFileName.c_str());  return 0; }

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
			
			getRelAbundance(lookup, out);
			
			for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }
		}
	
		//reset groups parameter
		globaldata->Groups.clear();  
		delete input; globaldata->ginput = NULL;
		delete read;
		out.close();
		
		if (m->control_pressed) { remove(outputFileName.c_str()); return 0;}
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		m->mothurOut(outputFileName); m->mothurOutEndLine();
		m->mothurOutEndLine();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "GetRelAbundCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************

int GetRelAbundCommand::getRelAbundance(vector<SharedRAbundVector*> thisLookUp, ofstream& out){
	try {
	
	
	}
	catch(exception& e) {
		m->errorOut(e, "GetRelAbundCommand", "getRelAbundance");
		exit(1);
	}
}
//**********************************************************************************************************************


