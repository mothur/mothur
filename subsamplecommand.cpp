/*
 *  subsamplecommand.cpp
 *  Mothur
 *
 *  Created by westcott on 10/27/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "subsamplecommand.h"
#include "sharedutilities.h"

//**********************************************************************************************************************
vector<string> SubSampleCommand::getValidParameters(){	
	try {
		string Array[] =  {"fasta", "group", "list","shared","rabund", "name","sabund","size","groups","label","outputdir","inputdir"};
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
		outputTypes["fasta"] = tempOutNames;
		outputTypes["name"] = tempOutNames;
		outputTypes["group"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "SubSampleCommand", "GetRelAbundCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> SubSampleCommand::getRequiredParameters(){	
	try {
		string Array[] =  {"fasta","list","shared","rabund", "sabund","or"};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
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
		vector<string> myArray;
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
			string Array[] =  {"fasta", "group", "list","shared","rabund", "sabund","name","size","groups","label","outputdir","inputdir"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
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
			outputTypes["list"] = tempOutNames;
			outputTypes["rabund"] = tempOutNames;
			outputTypes["sabund"] = tempOutNames;
			outputTypes["fasta"] = tempOutNames;
			outputTypes["name"] = tempOutNames;
			outputTypes["group"] = tempOutNames;
					
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";	}
			
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
				
				it = parameters.find("fasta");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["fasta"] = inputDir + it->second;		}
				}
				
				it = parameters.find("shared");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["shared"] = inputDir + it->second;		}
				}
				
				it = parameters.find("group");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["group"] = inputDir + it->second;		}
				}
				
				it = parameters.find("sabund");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["sabund"] = inputDir + it->second;		}
				}
				
				it = parameters.find("rabund");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["rabund"] = inputDir + it->second;		}
				}
				
				it = parameters.find("name");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["name"] = inputDir + it->second;		}
				}
			}
			
			//check for required parameters
			listfile = validParameter.validFile(parameters, "list", true);
			if (listfile == "not open") { listfile = ""; abort = true; }
			else if (listfile == "not found") { listfile = ""; }	
			
			sabundfile = validParameter.validFile(parameters, "sabund", true);
			if (sabundfile == "not open") { sabundfile = ""; abort = true; }	
			else if (sabundfile == "not found") { sabundfile = ""; }
			
			rabundfile = validParameter.validFile(parameters, "rabund", true);
			if (rabundfile == "not open") { rabundfile = ""; abort = true; }	
			else if (rabundfile == "not found") { rabundfile = ""; }
			
			fastafile = validParameter.validFile(parameters, "fasta", true);
			if (fastafile == "not open") { fastafile = ""; abort = true; }	
			else if (fastafile == "not found") { fastafile = ""; }
			
			sharedfile = validParameter.validFile(parameters, "shared", true);
			if (sharedfile == "not open") { sharedfile = ""; abort = true; }	
			else if (sharedfile == "not found") { sharedfile = ""; }
			
			namefile = validParameter.validFile(parameters, "name", true);
			if (namefile == "not open") { namefile = ""; abort = true; }	
			else if (namefile == "not found") { namefile = ""; }
			
			groupfile = validParameter.validFile(parameters, "group", true);
			if (groupfile == "not open") { groupfile = ""; abort = true; }	
			else if (groupfile == "not found") { groupfile = ""; }
			
			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  m->splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
			
			groups = validParameter.validFile(parameters, "groups", false);			
			if (groups == "not found") { groups = ""; pickedGroups = false; }
			else { 
				pickedGroups = true;
				m->splitAtDash(groups, Groups);
				globaldata->Groups = Groups;
			}
			
			string temp = validParameter.validFile(parameters, "size", false);		if (temp == "not found"){	temp = "0";		}
			convert(temp, size);  
			
			if ((namefile != "") && (fastafile == "")) { m->mothurOut("You may only use a namefile with a fastafile."); m->mothurOutEndLine(); abort = true; }
			
			if ((fastafile == "") && (listfile == "") && (sabundfile == "") && (rabundfile == "") && (sharedfile == "")) {
				m->mothurOut("You must provide a fasta, list, sabund, rabund or shared file as an input file."); m->mothurOutEndLine(); abort = true; }
			
			if (pickedGroups && ((groupfile == "") && (sharedfile == ""))) { 
				m->mothurOut("You cannot pick groups without a valid group file or shared file."); m->mothurOutEndLine(); abort = true; }
			
			if ((groupfile != "") && ((fastafile == "") && (listfile == ""))) { 
				m->mothurOut("Group file only valid with listfile or fastafile."); m->mothurOutEndLine(); abort = true; }
			
			if ((groupfile != "") && ((fastafile != "") && (listfile != ""))) { 
				m->mothurOut("A new group file can only be made from the subsample of a listfile or fastafile, not both. Please correct."); m->mothurOutEndLine(); abort = true; }
			
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
		m->mothurOut("The sub.sample command is designed to be used as a way to normalize your data, or create a smaller set from your original set.\n");
		m->mothurOut("The sub.sample command parameters are fasta, name, list, group, rabund, sabund, shared, groups, size and label.  You must provide a fasta, list, sabund, rabund or shared file as an input file.\n");
		m->mothurOut("The groups parameter allows you to specify which of the groups in your groupfile you would like included. The group names are separated by dashes.\n");
		m->mothurOut("The label parameter allows you to select what distance levels you would like, and are also separated by dashes.\n");
		m->mothurOut("The size parameter allows you indicate the size of your subsample.\n");
		m->mothurOut("The size parameter is not set: with shared file size=number of seqs in smallest sample, with all other files, 10% of number of seqs.\n");
		m->mothurOut("The sub.sample command should be in the following format: sub.sample(list=yourListFile, group=yourGroupFile, groups=yourGroups, label=yourLabels).\n");
		m->mothurOut("Example sub.sample(list=abrecovery.fn.list, group=abrecovery.groups, groups=B-C, size=20).\n");
		m->mothurOut("The default value for groups is all the groups in your groupfile, and all labels in your inputfile will be used.\n");
		m->mothurOut("The sub.sample command outputs a .subsample file.\n");
		m->mothurOut("Note: No spaces between parameter labels (i.e. groups), '=' and parameters (i.e.yourGroups).\n\n");

	}
	catch(exception& e) {
		m->errorOut(e, "SubSampleCommand", "help");
		exit(1);
	}
}

//**********************************************************************************************************************

SubSampleCommand::~SubSampleCommand(){}

//**********************************************************************************************************************

int SubSampleCommand::execute(){
	try {
	
		if (abort == true) { return 0; }
		
		if (sharedfile != "")	{   getSubSampleShared();	}
		if (m->control_pressed) {  for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); return 0; } }
		if (listfile != "")		{   getSubSampleList();		}
		if (m->control_pressed) {  for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); return 0; } }
		//if (rabund != "")		{   getSubSampleRabund();	}
		if (m->control_pressed) {  for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); return 0; } }
		//if (sabundfile != "")	{   getSubSampleSabund();	}
		if (m->control_pressed) {  for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); return 0; } }
		//if (fastafile != "")	{   getSubSampleFasta();	}
		if (m->control_pressed) {  for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); return 0; } }
			
				
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SubSampleCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
int SubSampleCommand::getSubSampleShared() {
	try {
		
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += m->hasPath(sharedfile);  }
		string outputFileName = thisOutputDir + m->getRootName(m->getSimpleName(sharedfile)) + "subsample" + m->getExtension(sharedfile);
		
		ofstream out;
		m->openOutputFile(outputFileName, out);
		outputTypes["shared"].push_back(outputFileName);  outputNames.push_back(outputFileName);
		
		InputData* input = new InputData(sharedfile, "sharedfile");
		vector<SharedRAbundVector*> lookup = input->getSharedRAbundVectors();
		string lastLabel = lookup[0]->getLabel();
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;
		
		
		//as long as you are not at the end of the file or done wih the lines you want
		while((lookup[0] != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			if (m->control_pressed) {  delete input; out.close(); return 0;  }
			
			if(allLines == 1 || labels.count(lookup[0]->getLabel()) == 1){			
				
				m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
				
				processShared(lookup, out);
				
				processedLabels.insert(lookup[0]->getLabel());
				userLabels.erase(lookup[0]->getLabel());
			}
			
			if ((m->anyLabelsToProcess(lookup[0]->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
				string saveLabel = lookup[0]->getLabel();
				
				for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }  
				
				lookup = input->getSharedRAbundVectors(lastLabel);
				m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
				
				processShared(lookup, out);
				
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
		
		
		if (m->control_pressed) {  out.close(); return 0;  }
		
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
			
			processShared(lookup, out);
			
			for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }
		}
		
		delete input;
		out.close();  
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "SubSampleCommand", "getSubSampleShared");
		exit(1);
	}
}
//**********************************************************************************************************************
int SubSampleCommand::processShared(vector<SharedRAbundVector*>& thislookup, ofstream& out) {
	try {
		
		//if (pickedGroups) { eliminateZeroOTUS(thislookup); }
		
		if (size == 0) { //user has not set size, set size = smallest samples size
			size = thislookup[0]->getNumSeqs();
			for (int i = 1; i < thislookup.size(); i++) {
				int thisSize = thislookup[i]->getNumSeqs();
				
				if (thisSize < size) {	size = thisSize;	}
			}
		}
		
		int numBins = thislookup[0]->getNumBins();
		for (int i = 0; i < thislookup.size(); i++) {		
			int thisSize = thislookup[i]->getNumSeqs();
			
			if (thisSize != size) {
				
				string thisgroup = thislookup[i]->getGroup();
				
				OrderVector* order = new OrderVector();
				for(int p=0;p<numBins;p++){
					for(int j=0;j<thislookup[i]->getAbundance(p);j++){
						order->push_back(p);
					}
				}
				random_shuffle(order->begin(), order->end());
				
				SharedRAbundVector* temp = new SharedRAbundVector(numBins);
				temp->setLabel(thislookup[i]->getLabel());
				temp->setGroup(thislookup[i]->getGroup());
				
				delete thislookup[i];
				thislookup[i] = temp;
				
				
				for (int j = 0; j < size; j++) {
					
					if (m->control_pressed) { delete order; return 0; }
					
					//get random number to sample from order between 0 and thisSize-1.
					int myrand = (int)((float)(rand()) / (RAND_MAX / (thisSize-1) + 1));
					
					int bin = order->get(myrand);
					
					int abund = thislookup[i]->getAbundance(bin);
					thislookup[i]->set(bin, (abund+1), thisgroup);
				}	
				delete order;
			}
		}
		
		//subsampling may have created some otus with no sequences in them
		eliminateZeroOTUS(thislookup);
		
		if (m->control_pressed) { return 0; }
		
		for (int i = 0; i < thislookup.size(); i++) {
			out << thislookup[i]->getLabel() << '\t' << thislookup[i]->getGroup() << '\t';
			thislookup[i]->print(out);
		}
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "SubSampleCommand", "processShared");
		exit(1);
	}
}			
//**********************************************************************************************************************
int SubSampleCommand::getSubSampleList() {
	try {
		
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += m->hasPath(listfile);  }
		string outputFileName = thisOutputDir + m->getRootName(m->getSimpleName(listfile)) + "subsample" + m->getExtension(listfile);
		
		ofstream out;
		m->openOutputFile(outputFileName, out);
		outputTypes["list"].push_back(outputFileName);  outputNames.push_back(outputFileName);
		
		InputData* input;
		//if you have a groupfile you want to read a shared list
		if (groupfile != "") {
			
			GroupMap* groupMap = new GroupMap(groupfile);
			groupMap->readMap();
			
			//takes care of user setting groupNames that are invalid or setting groups=all
			SharedUtil* util = new SharedUtil();
			util->setGroups(Groups, groupMap->namesOfGroups);
			delete util;
			
			//create outputfiles
			string groupOutputDir = outputDir;
			if (outputDir == "") {  groupOutputDir += m->hasPath(groupfile);  }
			string groupOutputFileName = groupOutputDir + m->getRootName(m->getSimpleName(groupfile)) + "subsample" + m->getExtension(groupfile);
			
			ofstream outGroup;
			m->openOutputFile(groupOutputFileName, outGroup);
			outputTypes["group"].push_back(groupOutputFileName);  outputNames.push_back(groupOutputFileName);
			
			globaldata->setGroupFile(groupfile); //shared list needs this
			
			input = new InputData(listfile, "list");
			ListVector* list = input->getListVector();
			string lastLabel = list->getLabel();
			
			//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
			set<string> processedLabels;
			set<string> userLabels = labels;
			
			//file mismatch quit
			if (list->getNumSeqs() != groupMap->getNumSeqs()) { 
				m->mothurOut("[ERROR]: your list file contains " + toString(list->getNumSeqs()) + " sequences, and your groupfile contains " + toString(groupMap->getNumSeqs()) + ", please correct."); 
				m->mothurOutEndLine();
				delete groupMap;
				delete list;
				delete input;
				out.close();
				outGroup.close();
				return 0;
			}
			
			//as long as you are not at the end of the file or done wih the lines you want
			while((list != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
				
				if (m->control_pressed) {  delete list; delete input; delete groupMap; out.close(); outGroup.close(); return 0;  }
				
				if(allLines == 1 || labels.count(list->getLabel()) == 1){			
					
					m->mothurOut(list->getLabel()); m->mothurOutEndLine();
					
					processListGroup(list, groupMap, out, outGroup);
					
					processedLabels.insert(list->getLabel());
					userLabels.erase(list->getLabel());
				}
				
				if ((m->anyLabelsToProcess(list->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
					string saveLabel = list->getLabel();
					
					delete list; 
					
					list = input->getListVector(lastLabel);
					m->mothurOut(list->getLabel()); m->mothurOutEndLine();
					
					processListGroup(list, groupMap, out, outGroup);
					
					processedLabels.insert(list->getLabel());
					userLabels.erase(list->getLabel());
					
					//restore real lastlabel to save below
					list->setLabel(saveLabel);
				}
				
				lastLabel = list->getLabel();
				
				delete list; list = NULL;
				
				//get next line to process
				list = input->getListVector();				
			}
			
			
			if (m->control_pressed) {  if (list != NULL) { delete list; } delete input; delete groupMap; out.close(); outGroup.close(); return 0;  }
			
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
				if (list != NULL) { delete list; }
				
				list = input->getListVector(lastLabel);
				
				m->mothurOut(list->getLabel()); m->mothurOutEndLine();
				
				processListGroup(list, groupMap, out, outGroup);
				
				delete list; list = NULL;
			}
			
			out.close();  outGroup.close();
			if (list != NULL) { delete list; }
			delete groupMap;
			
		}else {
			//need to complete
		}
		
		
		delete input;
						
		return 0;
 
	}
	catch(exception& e) {
		m->errorOut(e, "SubSampleCommand", "getSubSampleList");
		exit(1);
	}
}
//**********************************************************************************************************************
int SubSampleCommand::processListGroup(ListVector*& list, GroupMap*& groupMap, ofstream& out, ofstream& outGroup) {
	try {
				
		if (size == 0) { //user has not set size, set size = 10% samples size
			size = int (list->getNumSeqs() * 0.10);
		}
		
		int numBins = list->getNumBins();
		int thisSize = list->getNumSeqs();
		
		if (size > thisSize) { m->mothurOut("Your list file only contains " + toString(thisSize) + " sequences. Setting size to " + toString(thisSize) + "."); m->mothurOutEndLine();
			size = thisSize;
		}
		
		vector<nameToBin> seqs;
		for (int i = 0; i < numBins; i++) {
			string names = list->get(i);
			
			//parse names
			string individual = "";
			int length = names.length();
			for(int j=0;j<length;j++){
				if(names[j] == ','){
					nameToBin temp(individual, i);
					seqs.push_back(temp);
					individual = "";				
				}
				else{
					individual += names[j];
				}
			}
			nameToBin temp(individual, i);
			seqs.push_back(temp);
		}
		
					
		ListVector* temp = new ListVector(numBins);
		temp->setLabel(list->getLabel());
		
		delete list; 
		list = temp;
		
		set<int> alreadySelected; //dont want repeat sequence names added
		alreadySelected.insert(-1);
		for (int j = 0; j < size; j++) {
			
			if (m->control_pressed) { return 0; }
			
			//get random sequence to add, making sure we have not already added it
			int myrand = -1;
			while (alreadySelected.count(myrand) != 0) {
				myrand = (int)((float)(rand()) / (RAND_MAX / (thisSize-1) + 1));
			}
			alreadySelected.insert(myrand);
			
			//update new list
			string oldNames = temp->get(seqs[myrand].bin);
			if (oldNames == "") { oldNames += seqs[myrand].name; }
			else { oldNames += "," + seqs[myrand].name; }
			
			temp->set(seqs[myrand].bin, oldNames);
			
			//update group file
			string group = groupMap->getGroup(seqs[myrand].name);
			if (group == "not found") { m->mothurOut("[ERROR]: " + seqs[myrand].name + " is not in your groupfile. please correct."); m->mothurOutEndLine(); group = "NOTFOUND"; }
			
			outGroup << seqs[myrand].name << '\t' << group << endl;
		}	

		if (m->control_pressed) { return 0; }
		
		list->print(out);
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "SubSampleCommand", "processListGroup");
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
		thislookup.clear();
		
		thislookup = newLookup;
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "SubSampleCommand", "eliminateZeroOTUS");
		exit(1);
	}
}

//**********************************************************************************************************************



