/*
 *  subsamplecommand.cpp
 *  Mothur
 *
 *  Created by westcott on 10/27/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "subsamplecommand.h"

//**********************************************************************************************************************
vector<string> SubSampleCommand::getValidParameters(){	
	try {
		string Array[] =  {"groups","label","outputdir","inputdir"};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "SubSampleCommand", "getValidParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
SubSampleCommand::SubSampleCommand(){	
	try {
		abort = true;
		//initialize outputTypes
		vector<string> tempOutNames;
		outputTypes["shared"] = tempOutNames;
		outputTypes["list"] = tempOutNames;
		outputTypes["rabund"] = tempOutNames;
		outputTypes["sabund"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "SubSampleCommand", "GetRelAbundCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> SubSampleCommand::getRequiredParameters(){	
	try {
		vector<string> myArray;
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "SubSampleCommand", "getRequiredParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> SubSampleCommand::getRequiredFiles(){	
	try {
		string Array[] =  {"shared","list","rabund","sabund","or"};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "SubSampleCommand", "getRequiredFiles");
		exit(1);
	}
}
//**********************************************************************************************************************
SubSampleCommand::SubSampleCommand(string option) {
	try {
		globaldata = GlobalData::getInstance();
		abort = false;
		allLines = 1;
		labels.clear();
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string AlignArray[] =  {"groups","label","outputdir","inputdir"};
			vector<string> myArray (AlignArray, AlignArray+(sizeof(AlignArray)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			
			//check to make sure all parameters are valid for command
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["shared"] = tempOutNames;
			outputTypes["list"] = tempOutNames;
			outputTypes["rabund"] = tempOutNames;
			outputTypes["sabund"] = tempOutNames;
					
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += m->hasPath(globaldata->inputFileName); //if user entered a file with a path then preserve it	
			}
			
			//make sure the user has already run the read.otu command
			if ((globaldata->getSharedFile() == "") && (globaldata->getListFile() == "") && (globaldata->getRabundFile() == "") && (globaldata->getSabundFile() == "")) { m->mothurOut("You must read a list, sabund, rabund or shared file before you can use the sub.sample command."); m->mothurOutEndLine(); abort = true; }

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
			if (groups == "not found") { groups = ""; pickedGroups = false; }
			else { 
				pickedGroups = true;
				m->splitAtDash(groups, Groups);
				globaldata->Groups = Groups;
			}
			
		}

	}
	catch(exception& e) {
		m->errorOut(e, "SubSampleCommand", "SubSampleCommand");
		exit(1);
	}
}

//**********************************************************************************************************************

void SubSampleCommand::help(){
	try {
		m->mothurOut("The get.relabund command can only be executed after a successful read.otu command of a list and group or shared file.\n");
		m->mothurOut("The get.relabund command parameters are groups, scale and label.  No parameters are required.\n");
		m->mothurOut("The groups parameter allows you to specify which of the groups in your groupfile you would like included. The group names are separated by dashes.\n");
		m->mothurOut("The label parameter allows you to select what distance levels you would like, and are also separated by dashes.\n");
		m->mothurOut("The scale parameter allows you to select what scale you would like to use. Choices are totalgroup, totalotu, averagegroup, averageotu, default is totalgroup.\n");
		m->mothurOut("The get.relabund command should be in the following format: get.relabund(groups=yourGroups, label=yourLabels).\n");
		m->mothurOut("Example get.relabund(groups=A-B-C, scale=averagegroup).\n");
		m->mothurOut("The default value for groups is all the groups in your groupfile, and all labels in your inputfile will be used.\n");
		m->mothurOut("The get.relabund command outputs a .relabund file.\n");
		m->mothurOut("Note: No spaces between parameter labels (i.e. groups), '=' and parameters (i.e.yourGroups).\n\n");

	}
	catch(exception& e) {
		m->errorOut(e, "SubSampleCommand", "help");
		exit(1);
	}
}

//**********************************************************************************************************************

SubSampleCommand::~SubSampleCommand(){
}

//**********************************************************************************************************************

int SubSampleCommand::execute(){
	try {
	
		if (abort == true) { return 0; }
		
		string outputFileName = outputDir + m->getRootName(m->getSimpleName(globaldata->inputFileName)) + "subsample" +  m->getExtension(globaldata->inputFileName);
		ofstream out;
		m->openOutputFile(outputFileName, out);
		out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);
		
		string format = globaldata->getFormat();
		
		read = new ReadOTUFile(globaldata->inputFileName);	
		read->read(&*globaldata); 
		input = globaldata->ginput;
		
		if (format == "sharedfile") {
			lookup = input->getSharedRAbundVectors();
			outputTypes["shared"].push_back(outputFileName);
			getSubSampleShared(lookup, out);
		}else if (format == "list") { 
			list = globaldata->glist;
			outputTypes["list"].push_back(outputFileName);
			//getSubSamplesList();
		}else if (format == "rabund") { 
			rabund = globaldata->rabund;
			outputTypes["rabund"].push_back(outputFileName);
			//getSubSamplesRabund();
		
		}else if (format == "sabund") { 
			sabund = globaldata->sabund;
			outputTypes["sabund"].push_back(outputFileName);
			//getSubSamplesSabund();
		}
		
		out.close();
					
		//reset groups parameter
		delete input; globaldata->ginput = NULL;
		delete read;
		
		if (m->control_pressed) { outputTypes.clear(); remove(outputFileName.c_str()); return 0;}
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		m->mothurOut(outputFileName); m->mothurOutEndLine(); outputNames.push_back(outputFileName); 
		m->mothurOutEndLine();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SubSampleCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
int SubSampleCommand::getSubSampleShared(vector<SharedRAbundVector*>& thislookup, ofstream& filename) {
	try {
	
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;

		string lastLabel = lookup[0]->getLabel();
	
		//as long as you are not at the end of the file or done wih the lines you want
		while((lookup[0] != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			if (m->control_pressed) {  return 0;  }
	
			if(allLines == 1 || labels.count(lookup[0]->getLabel()) == 1){			

				m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
				
				//process lookup
				
								
																
				processedLabels.insert(lookup[0]->getLabel());
				userLabels.erase(lookup[0]->getLabel());
			}
			
			if ((m->anyLabelsToProcess(lookup[0]->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
				string saveLabel = lookup[0]->getLabel();
		
				for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }  
		
				lookup = input->getSharedRAbundVectors(lastLabel);
				m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
				
				//process lookup				
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
		
		
		if (m->control_pressed) {  return 0;  }

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
			
			//process lookup
			
			
			
			for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }
		}
	
		//reset groups parameter
		globaldata->Groups.clear();  
				
		return 0;
 
	}
	catch(exception& e) {
		m->errorOut(e, "SubSampleCommand", "getSubSampleShared");
		exit(1);
	}
}


//**********************************************************************************************************************
int SubSampleCommand::eliminateZeroOTUS(vector<SharedRAbundVector*>& thislookup) {
	try {
		
		vector<SharedRAbundVector*> newLookup;
		for (int i = 0; i < thislookup.size(); i++) {
			SharedRAbundVector* temp = new SharedRAbundVector();
			temp->setLabel(thislookup[i]->getLabel());
			temp->setGroup(thislookup[i]->getGroup());
			newLookup.push_back(temp);
		}
		
		//for each bin
		for (int i = 0; i < thislookup[0]->getNumBins(); i++) {
			if (m->control_pressed) { for (int j = 0; j < newLookup.size(); j++) {  delete newLookup[j];  } return 0; }
		
			//look at each sharedRabund and make sure they are not all zero
			bool allZero = true;
			for (int j = 0; j < thislookup.size(); j++) {
				if (thislookup[j]->getAbundance(i) != 0) { allZero = false;  break;  }
			}
			
			//if they are not all zero add this bin
			if (!allZero) {
				for (int j = 0; j < thislookup.size(); j++) {
					newLookup[j]->push_back(thislookup[j]->getAbundance(i), thislookup[j]->getGroup());
				}
			}
		}

		for (int j = 0; j < thislookup.size(); j++) {  delete thislookup[j];  }

		thislookup = newLookup;
		
		return 0;
 
	}
	catch(exception& e) {
		m->errorOut(e, "SubSampleCommand", "eliminateZeroOTUS");
		exit(1);
	}
}

//**********************************************************************************************************************


