/*
 *  getlistcountcommand.cpp
 *  Mothur
 *
 *  Created by westcott on 10/12/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */

#include "getlistcountcommand.h"

//**********************************************************************************************************************
vector<string> GetListCountCommand::getValidParameters(){	
	try {
		string Array[] =  {"list","label","sort","outputdir","inputdir"};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "GetListCountCommand", "getValidParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
GetListCountCommand::GetListCountCommand(){	
	try {
		abort = true;
		//initialize outputTypes
		vector<string> tempOutNames;
		outputTypes["otu"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "GetListCountCommand", "GetListCountCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> GetListCountCommand::getRequiredParameters(){	
	try {
		string Array[] =  {"list"};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "GetListCountCommand", "getRequiredParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> GetListCountCommand::getRequiredFiles(){	
	try {
		vector<string> myArray;
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "GetListCountCommand", "getRequiredFiles");
		exit(1);
	}
}
//**********************************************************************************************************************
GetListCountCommand::GetListCountCommand(string option)  {
	try {
		globaldata = GlobalData::getInstance();
		abort = false;
		allLines = 1;
		labels.clear();
				
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string AlignArray[] =  {"list","label","sort","outputdir","inputdir"};
			vector<string> myArray (AlignArray, AlignArray+(sizeof(AlignArray)/sizeof(string)));
			
			OptionParser parser(option);
			map<string, string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			map<string, string>::iterator it;
			
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["otu"] = tempOutNames;
		
			string ranRead = globaldata->getListFile();
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("list");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["list"] = inputDir + it->second;		}
				}
			}

			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";		}
			
			//check for required parameters
			listfile = validParameter.validFile(parameters, "list", true);
			if ((listfile == "not found") && (globaldata->getListFile() == ""))  { m->mothurOut("You must read a listfile before running the get.listcount command.");  m->mothurOutEndLine(); abort = true; }
			else if ((listfile == "not found") && (globaldata->getListFile() != "")) { listfile = globaldata->getListFile(); }
			else if (listfile == "not open") { abort = true; }	
			else { globaldata->setListFile(listfile); }
		
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			sort = validParameter.validFile(parameters, "sort", false);	  if (sort == "not found") { sort = "otu"; }
			if ((sort != "otu") && (sort != "name")) { m->mothurOut( sort + " is not a valid sort option. Options are otu and name. I will use otu."); m->mothurOutEndLine(); sort = "otu"; }
			
			label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  m->splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
			
			//if the user has not specified any labels use the ones from read.otu
			if ((label == "") && (ranRead != "")) {  
				allLines = globaldata->allLines; 
				labels = globaldata->labels; 
			}
		}
	}
	catch(exception& e) {
		m->errorOut(e, "GetListCountCommand", "GetListCountCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

void GetListCountCommand::help(){
	try {
		m->mothurOut("The get.otulist command can only be executed after a successful read.otu command of a listfile or providing a list file using the list parameter.\n");
		m->mothurOut("The get.otulist command parameters are list, sort and label.  No parameters are required.\n");
		m->mothurOut("The label parameter allows you to select what distance levels you would like a output files created for, and are separated by dashes.\n");
		m->mothurOut("The sort parameter allows you to select how you want the output displayed. Options are otu and name.\n");
		m->mothurOut("If otu is selected the output will be otu number followed by the list of names in that otu.\n");
		m->mothurOut("If name is selected the output will be a sequence name followed by its otu number.\n");
		m->mothurOut("The get.otulist command should be in the following format: get.otulist(list=yourlistFile, label=yourLabels).\n");
		m->mothurOut("Example get.otulist(list=amazon.fn.list, label=0.10).\n");
		m->mothurOut("The default value for label is all lines in your inputfile.\n");
		m->mothurOut("The get.otulist command outputs a .otu file for each distance you specify listing the bin number and the names of the sequences in that bin.\n");
		m->mothurOut("Note: No spaces between parameter labels (i.e. list), '=' and parameters (i.e.yourListFile).\n\n");
	}
	catch(exception& e) {
		m->errorOut(e, "GetListCountCommand", "help");
		exit(1);
	}
}

//**********************************************************************************************************************

GetListCountCommand::~GetListCountCommand(){}

//**********************************************************************************************************************

int GetListCountCommand::execute(){
	try {
		if (abort == true) {	return 0;	}

		globaldata->setFormat("list");
		
		//read list file
		read = new ReadOTUFile(listfile);	
		read->read(&*globaldata); 
		
		input = globaldata->ginput;
		list = globaldata->gListVector;
		string lastLabel = list->getLabel();

		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;
		
		if (m->control_pressed) { 
			delete read;
			delete input;
			delete list;
			globaldata->gListVector = NULL;  
			for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());	}
			return 0; 
		}
		
		while((list != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			
			if(allLines == 1 || labels.count(list->getLabel()) == 1){
			
				process(list);
				
				if (m->control_pressed) { 
					delete read;
					delete input;
					delete list;
					globaldata->gListVector = NULL;  
					for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());	} outputTypes.clear();
					return 0; 
				}
							
				processedLabels.insert(list->getLabel());
				userLabels.erase(list->getLabel());
			}
			
			if ((m->anyLabelsToProcess(list->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
				string saveLabel = list->getLabel();
				
				delete list;
				list = input->getListVector(lastLabel);
				
				process(list);
				
				if (m->control_pressed) { 
					delete read;
					delete input;
					delete list;
					globaldata->gListVector = NULL;  
					for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());	} outputTypes.clear();
					return 0; 
				}
													
				processedLabels.insert(list->getLabel());
				userLabels.erase(list->getLabel());
				
				//restore real lastlabel to save below
				list->setLabel(saveLabel);
			}
			
			lastLabel = list->getLabel();			
			
			delete list;
			list = input->getListVector();
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
				m->mothurOut(". Please refer to " + lastLabel + ".");  m->mothurOutEndLine();
			}
		}
		
		//run last label if you need to
		if (needToRun == true)  {
			if (list != NULL) {		delete list;	}
			list = input->getListVector(lastLabel);
				
			process(list);	
			
			if (m->control_pressed) { 
					delete read;
					delete input;
					delete list;
					globaldata->gListVector = NULL;  
					for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());	} outputTypes.clear();
					return 0; 
			}
			
			delete list;  
		}
		
		delete read;
		delete input;
		globaldata->gListVector = NULL;
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "GetListCountCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************
//return 1 if error, 0 otherwise
void GetListCountCommand::process(ListVector* list) {
	try {
		string binnames;
		if (outputDir == "") { outputDir += m->hasPath(listfile); }
		string outputFileName = outputDir + m->getRootName(m->getSimpleName(listfile)) + list->getLabel() + ".otu";
		m->openOutputFile(outputFileName, out);
		outputNames.push_back(outputFileName); outputTypes["otu"].push_back(outputFileName);
		
		m->mothurOut(list->getLabel()); m->mothurOutEndLine();
		
		//for each bin in the list vector
		for (int i = 0; i < list->getNumBins(); i++) {
			if (m->control_pressed) { break; }
			
			binnames = list->get(i);
			
			if (sort == "otu") {
				out << i+1 << '\t' << binnames << endl;
			}else{ //sort = name
				vector<string> names;
				m->splitAtComma(binnames, names);
				
				for (int j = 0; j < names.size(); j++) {
					out << names[j] << '\t' << i+1 << endl;
				}
			}
		}
		
		out.close();
	}
	catch(exception& e) {
		m->errorOut(e, "GetListCountCommand", "process");
		exit(1);
	}
}
//**********************************************************************************************************************


