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
vector<string> SubSampleCommand::setParameters(){	
	try {		
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "FLSSR", "none",false,false); parameters.push_back(pfasta);
		CommandParameter pname("name", "InputTypes", "", "", "none", "none", "none",false,false); parameters.push_back(pname);
		CommandParameter pgroup("group", "InputTypes", "", "", "none", "none", "none",false,false); parameters.push_back(pgroup);
		CommandParameter plist("list", "InputTypes", "", "", "none", "FLSSR", "none",false,false); parameters.push_back(plist);
		CommandParameter pshared("shared", "InputTypes", "", "", "none", "FLSSR", "none",false,false); parameters.push_back(pshared);
		CommandParameter prabund("rabund", "InputTypes", "", "", "none", "FLSSR", "none",false,false); parameters.push_back(prabund);
		CommandParameter psabund("sabund", "InputTypes", "", "", "none", "FLSSR", "none",false,false); parameters.push_back(psabund);
		CommandParameter plabel("label", "String", "", "", "", "", "",false,false); parameters.push_back(plabel);
		CommandParameter pgroups("groups", "String", "", "", "", "", "",false,false); parameters.push_back(pgroups);
		CommandParameter psize("size", "Number", "", "0", "", "", "",false,false); parameters.push_back(psize);
		CommandParameter ppersample("persample", "Boolean", "", "F", "", "", "",false,false); parameters.push_back(ppersample);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "SubSampleCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string SubSampleCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The sub.sample command is designed to be used as a way to normalize your data, or create a smaller set from your original set.\n";
		helpString += "The sub.sample command parameters are fasta, name, list, group, rabund, sabund, shared, groups, size, persample and label.  You must provide a fasta, list, sabund, rabund or shared file as an input file.\n";
		helpString += "The namefile is only used with the fasta file, not with the listfile, because the list file should contain all sequences.\n";
		helpString += "The groups parameter allows you to specify which of the groups in your groupfile you would like included. The group names are separated by dashes.\n";
		helpString += "The label parameter allows you to select what distance levels you would like, and are also separated by dashes.\n";
		helpString += "The size parameter allows you indicate the size of your subsample.\n";
		helpString += "The persample parameter allows you indicate you want to select subsample of the same size from each of your groups, default=false. It is only used with the list and fasta files if a groupfile is given.\n";
		helpString += "persample=false will select a random set of sequences of the size you select, but the number of seqs from each group may differ.\n";
		helpString += "The size parameter is not set: with shared file size=number of seqs in smallest sample, with all other files if a groupfile is given and persample=true, then size=number of seqs in smallest sample, otherwise size=10% of number of seqs.\n";
		helpString += "The sub.sample command should be in the following format: sub.sample(list=yourListFile, group=yourGroupFile, groups=yourGroups, label=yourLabels).\n";
		helpString += "Example sub.sample(list=abrecovery.fn.list, group=abrecovery.groups, groups=B-C, size=20).\n";
		helpString += "The default value for groups is all the groups in your groupfile, and all labels in your inputfile will be used.\n";
		helpString += "The sub.sample command outputs a .subsample file.\n";
		helpString += "Note: No spaces between parameter labels (i.e. groups), '=' and parameters (i.e.yourGroups).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "SubSampleCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
SubSampleCommand::SubSampleCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
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
SubSampleCommand::SubSampleCommand(string option) {
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
				m->Groups = Groups;
			}
			
			string temp = validParameter.validFile(parameters, "size", false);		if (temp == "not found"){	temp = "0";		}
			convert(temp, size);  
			
			temp = validParameter.validFile(parameters, "persample", false);		if (temp == "not found"){	temp = "f";		}
			persample = m->isTrue(temp);
			
			if (groupfile == "") { persample = false; }
			
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

int SubSampleCommand::execute(){
	try {
	
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		if (sharedfile != "")	{   getSubSampleShared();	}
		if (m->control_pressed) {  for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); return 0; } }
		
		if (listfile != "")		{   getSubSampleList();		}
		if (m->control_pressed) {  for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); return 0; } }
		
		if (rabundfile != "")	{   getSubSampleRabund();	}
		if (m->control_pressed) {  for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); return 0; } }
		
		if (sabundfile != "")	{   getSubSampleSabund();	}
		if (m->control_pressed) {  for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); return 0; } }
		
		if (fastafile != "")	{   getSubSampleFasta();	}
		if (m->control_pressed) {  for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); return 0; } }
			
		//set fasta file as new current fastafile
		string current = "";
		itTypes = outputTypes.find("fasta");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setFastaFile(current); }
		}
		
		itTypes = outputTypes.find("name");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setNameFile(current); }
		}
		
		itTypes = outputTypes.find("group");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setGroupFile(current); }
		}
		
		itTypes = outputTypes.find("list");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setListFile(current); }
		}
		
		itTypes = outputTypes.find("shared");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setSharedFile(current); }
		}
		
		itTypes = outputTypes.find("rabund");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setRabundFile(current); }
		}
		
		itTypes = outputTypes.find("sabund");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setSabundFile(current); }
		}
		
		
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
int SubSampleCommand::getSubSampleFasta() {
	try {
		
		if (namefile != "") { readNames(); }	//fills names with all names in namefile.
		else { getNames(); }//no name file, so get list of names to pick from
		
		GroupMap* groupMap;
		if (groupfile != "") {
			
			groupMap = new GroupMap(groupfile);
			groupMap->readMap();
			
			//takes care of user setting groupNames that are invalid or setting groups=all
			SharedUtil* util = new SharedUtil();
			util->setGroups(Groups, groupMap->namesOfGroups);
			delete util;
			
			//file mismatch quit
			if (names.size() != groupMap->getNumSeqs()) { 
				m->mothurOut("[ERROR]: your fasta file contains " + toString(names.size()) + " sequences, and your groupfile contains " + toString(groupMap->getNumSeqs()) + ", please correct."); 
				m->mothurOutEndLine();
				delete groupMap;
				return 0;
			}			
		}	
		
		if (m->control_pressed) { return 0; }
		
		
		//make sure that if your picked groups size is not too big
		int thisSize = names.size();
		if (persample) { 
			if (size == 0) { //user has not set size, set size = smallest samples size
				size = groupMap->getNumSeqs(Groups[0]);
				for (int i = 1; i < Groups.size(); i++) {
					int thisSize = groupMap->getNumSeqs(Groups[i]);
					
					if (thisSize < size) {	size = thisSize;	}
				}
			}else { //make sure size is not too large
				vector<string> newGroups;
				for (int i = 0; i < Groups.size(); i++) {
					int thisSize = groupMap->getNumSeqs(Groups[i]);
					
					if (thisSize >= size) {	newGroups.push_back(Groups[i]);	}
					else {  m->mothurOut("You have selected a size that is larger than " + Groups[i] + " number of sequences, removing " + Groups[i] + "."); m->mothurOutEndLine(); }
				}
				Groups = newGroups;
			}
			
			m->mothurOut("Sampling " + toString(size) + " from each group."); m->mothurOutEndLine();			
		}else {
			if (pickedGroups) {
				int total = 0;
				for(int i = 0; i < Groups.size(); i++) {
					total += groupMap->getNumSeqs(Groups[i]);
				}
				
				if (size == 0) { //user has not set size, set size = 10% samples size
					size = int (total * 0.10);
				}
				
				if (total < size) { 
					if (size != 0) { 
						m->mothurOut("Your size is too large for the number of groups you selected. Adjusting to " + toString(int (total * 0.10)) + "."); m->mothurOutEndLine();
					}
					size = int (total * 0.10);
				}
				
				m->mothurOut("Sampling " + toString(size) + " from " + toString(total) + "."); m->mothurOutEndLine();
			}
			
			if (size == 0) { //user has not set size, set size = 10% samples size
				size = int (names.size() * 0.10);
			}
			
			if (size > thisSize) { m->mothurOut("Your fasta file only contains " + toString(thisSize) + " sequences. Setting size to " + toString(thisSize) + "."); m->mothurOutEndLine();
				size = thisSize;
			}
			
			if (!pickedGroups) { m->mothurOut("Sampling " + toString(size) + " from " + toString(thisSize) + "."); m->mothurOutEndLine(); }

		}
		random_shuffle(names.begin(), names.end());
		
		set<string> subset; //dont want repeat sequence names added
		if (persample) {
			for (int i = 0; i < Groups.size(); i++) {
				
				//randomly select a subset of those names from this group to include in the subsample
				for (int j = 0; j < size; j++) {
					
					if (m->control_pressed) { return 0; }
					
					//get random sequence to add, making sure we have not already added it
					bool done = false;
					int myrand;
					while (!done) {
						myrand = int((float)(thisSize) * (float)(rand()) / ((float)RAND_MAX+1.0));
						
						if (subset.count(names[myrand]) == 0)  { 
							
							string group = groupMap->getGroup(names[myrand]);
							if (group == "not found") { m->mothurOut("[ERROR]: " + names[myrand] + " is not in your groupfile. please correct."); m->mothurOutEndLine(); group = "NOTFOUND"; }
								
							if (group == Groups[i]) {	subset.insert(names[myrand]); break;	}
						}
					}
				}
			}
		}else {
			
			//randomly select a subset of those names to include in the subsample
			for (int j = 0; j < size; j++) {
				
				if (m->control_pressed) { return 0; }
				
				//get random sequence to add, making sure we have not already added it
				bool done = false;
				int myrand;
				while (!done) {
					myrand = int((float)(thisSize) * (float)(rand()) / ((float)RAND_MAX+1.0));
					
					if (subset.count(names[myrand]) == 0)  { 
						
						if (groupfile != "") { //if there is a groupfile given fill in group info
							string group = groupMap->getGroup(names[myrand]);
							if (group == "not found") { m->mothurOut("[ERROR]: " + names[myrand] + " is not in your groupfile. please correct."); m->mothurOutEndLine(); group = "NOTFOUND"; }
							
							if (pickedGroups) { //if hte user picked groups, we only want to keep the names of sequences from those groups
								if (m->inUsersGroups(group, Groups)) {
									subset.insert(names[myrand]); break;
								}
							}else{
								subset.insert(names[myrand]); break;
							}
						}else{ //save everyone, group
							subset.insert(names[myrand]); break;
						}					
					}
				}
			}	
		}
		
		if (subset.size() == 0) {  m->mothurOut("The size you selected is too large, skipping fasta file."); m->mothurOutEndLine();  return 0; }
		
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += m->hasPath(fastafile);  }
		string outputFileName = thisOutputDir + m->getRootName(m->getSimpleName(fastafile)) + "subsample" + m->getExtension(fastafile);
		
		ofstream out;
		m->openOutputFile(outputFileName, out);
		outputTypes["fasta"].push_back(outputFileName);  outputNames.push_back(outputFileName);
		
		//read through fasta file outputting only the names on the subsample list
		ifstream in;
		m->openInputFile(fastafile, in);
		
		string thisname;
		int count = 0;
		map<string, vector<string> >::iterator itNameMap;
		
		while(!in.eof()){
			
			if (m->control_pressed) { in.close(); out.close();  return 0; }
			
			Sequence currSeq(in);
			thisname = currSeq.getName();
			
			if (thisname != "") {
				
				//does the subset contain a sequence that this sequence represents
				itNameMap = nameMap.find(thisname);
				if (itNameMap != nameMap.end()) {
					vector<string> nameRepresents = itNameMap->second;
				
					for (int i = 0; i < nameRepresents.size(); i++){
						if (subset.count(nameRepresents[i]) != 0) {
							out << ">" << nameRepresents[i] << endl << currSeq.getAligned() << endl;
							count++;
						}
					}
				}else{
					m->mothurOut("[ERROR]: " + thisname + " is not in your namefile, please correct."); m->mothurOutEndLine();
				}
			}
			m->gobble(in);
		}
		in.close();	
		out.close();
		
		if (count != subset.size()) {
			m->mothurOut("[ERROR]: The subset selected contained " + toString(subset.size()) + " sequences, but I only found " + toString(count) + " of those in the fastafile."); m->mothurOutEndLine();
		}
		
		//if a groupfile is provided read through the group file only outputting the names on the subsample list
		if (groupfile != "") {
			
			string groupOutputDir = outputDir;
			if (outputDir == "") {  groupOutputDir += m->hasPath(groupfile);  }
			string groupOutputFileName = groupOutputDir + m->getRootName(m->getSimpleName(groupfile)) + "subsample" + m->getExtension(groupfile);
			
			ofstream outGroup;
			m->openOutputFile(groupOutputFileName, outGroup);
			outputTypes["group"].push_back(groupOutputFileName);  outputNames.push_back(groupOutputFileName);
			
			ifstream inGroup;
			m->openInputFile(groupfile, inGroup);
			string name, group;
			
			while(!inGroup.eof()){
				
				if (m->control_pressed) { inGroup.close(); outGroup.close(); return 0; }
				
				inGroup >> name;	m->gobble(inGroup);			//read from first column
				inGroup >> group;			//read from second column
				
				//if this name is in the accnos file
				if (subset.count(name) != 0) {
					outGroup << name << '\t' << group << endl;
					subset.erase(name);
				}
				
				m->gobble(inGroup);
			}
			inGroup.close();
			outGroup.close();	
			
			//sanity check
			if (subset.size() != 0) {  
				m->mothurOut("Your groupfile does not match your fasta file."); m->mothurOutEndLine();
				for (set<string>::iterator it = subset.begin(); it != subset.end(); it++) {
					m->mothurOut("[ERROR]: " + *it + " is missing from your groupfile."); m->mothurOutEndLine();
				}
			}
		}
			
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "SubSampleCommand", "getSubSampleFasta");
		exit(1);
	}
}
//**********************************************************************************************************************
int SubSampleCommand::getNames() {
	try {
		
		ifstream in;
		m->openInputFile(fastafile, in);
		
		string thisname;
		while(!in.eof()){
			
			if (m->control_pressed) { in.close(); return 0; }
			
			Sequence currSeq(in);
			thisname = currSeq.getName();
			
			if (thisname != "") {
				vector<string> temp; temp.push_back(thisname);
				nameMap[thisname] = temp;
				names.push_back(thisname);
			}
			m->gobble(in);
		}
		in.close();	
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "SubSampleCommand", "getNames");
		exit(1);
	}
}	
//**********************************************************************************************************************
int SubSampleCommand::readNames() {
	try {
		
		ifstream in;
		m->openInputFile(namefile, in);
		
		string thisname, repnames;
		map<string, vector<string> >::iterator it;
		
		while(!in.eof()){
			
			if (m->control_pressed) { in.close(); return 0; }
			
			in >> thisname;		m->gobble(in);		//read from first column
			in >> repnames;			//read from second column
			
			it = nameMap.find(thisname);
			if (it == nameMap.end()) {
				
				vector<string> splitRepNames;
				m->splitAtComma(repnames, splitRepNames);
				
				nameMap[thisname] = splitRepNames;	
				for (int i = 0; i < splitRepNames.size(); i++) { names.push_back(splitRepNames[i]); }
				
			}else{	m->mothurOut(thisname + " is already in namesfile. I will use first definition."); m->mothurOutEndLine();  }
			
			m->gobble(in);
		}
		in.close();	
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "SubSampleCommand", "readNames");
		exit(1);
	}
}		
//**********************************************************************************************************************
int SubSampleCommand::getSubSampleShared() {
	try {
		
		InputData* input = new InputData(sharedfile, "sharedfile");
		vector<SharedRAbundVector*> lookup = input->getSharedRAbundVectors();
		string lastLabel = lookup[0]->getLabel();
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;
		
		if (size == 0) { //user has not set size, set size = smallest samples size
			size = lookup[0]->getNumSeqs();
			for (int i = 1; i < lookup.size(); i++) {
				int thisSize = lookup[i]->getNumSeqs();
				
				if (thisSize < size) {	size = thisSize;	}
			}
		}else {
			m->Groups.clear();
			vector<SharedRAbundVector*> temp;
			for (int i = 0; i < lookup.size(); i++) {
				if (lookup[i]->getNumSeqs() < size) { 
					m->mothurOut(lookup[i]->getGroup() + " contains " + toString(lookup[i]->getNumSeqs()) + ". Eliminating."); m->mothurOutEndLine();
					delete lookup[i];
				}else { 
					m->Groups.push_back(lookup[i]->getGroup()); 
					temp.push_back(lookup[i]);
				}
			} 
			lookup = temp;
			Groups = m->Groups;
		}
		
		if (lookup.size() == 0) {  m->mothurOut("The size you selected is too large, skipping shared file."); m->mothurOutEndLine(); delete input; return 0; }
		
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += m->hasPath(sharedfile);  }
		string outputFileName = thisOutputDir + m->getRootName(m->getSimpleName(sharedfile)) + "subsample" + m->getExtension(sharedfile);
		
		ofstream out;
		m->openOutputFile(outputFileName, out);
		outputTypes["shared"].push_back(outputFileName);  outputNames.push_back(outputFileName);
		
		
		m->mothurOut("Sampling " + toString(size) + " from each group."); m->mothurOutEndLine();
		
		//as long as you are not at the end of the file or done wih the lines you want
		while((lookup[0] != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			if (m->control_pressed) {  delete input; for (int i = 0; i < lookup.size(); i++) {  delete lookup[i]; lookup[i] = NULL; } out.close(); return 0;  }
			
			if(allLines == 1 || labels.count(lookup[0]->getLabel()) == 1){			
				
				m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
				
				if (!m->printedHeaders) { lookup[0]->printHeaders(out); }
				processShared(lookup, out);
				
				processedLabels.insert(lookup[0]->getLabel());
				userLabels.erase(lookup[0]->getLabel());
			}
			
			if ((m->anyLabelsToProcess(lookup[0]->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
				string saveLabel = lookup[0]->getLabel();
				
				for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }  
				
				lookup = input->getSharedRAbundVectors(lastLabel);
				m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
				
				if (!m->printedHeaders) { lookup[0]->printHeaders(out); }
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
			
			if (!m->printedHeaders) { lookup[0]->printHeaders(out); }
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
					int myrand = int((float)(thisSize) * (float)(rand()) / ((float)RAND_MAX+1.0));
					
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
		
		InputData* input = new InputData(listfile, "list");
		ListVector* list = input->getListVector();
		string lastLabel = list->getLabel();
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;
		
		ofstream outGroup;
		GroupMap* groupMap;
		if (groupfile != "") {
			
			groupMap = new GroupMap(groupfile);
			groupMap->readMap();
			
			//takes care of user setting groupNames that are invalid or setting groups=all
			SharedUtil* util = new SharedUtil();
			util->setGroups(Groups, groupMap->namesOfGroups);
			delete util;
			
			//create outputfiles
			string groupOutputDir = outputDir;
			if (outputDir == "") {  groupOutputDir += m->hasPath(groupfile);  }
			string groupOutputFileName = groupOutputDir + m->getRootName(m->getSimpleName(groupfile)) + "subsample" + m->getExtension(groupfile);
			
			m->openOutputFile(groupOutputFileName, outGroup);
			outputTypes["group"].push_back(groupOutputFileName);  outputNames.push_back(groupOutputFileName);
			
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
		}
		
		//make sure that if your picked groups size is not too big
		if (persample) {
			if (size == 0) { //user has not set size, set size = smallest samples size
				size = groupMap->getNumSeqs(Groups[0]);
				for (int i = 1; i < Groups.size(); i++) {
					int thisSize = groupMap->getNumSeqs(Groups[i]);
					
					if (thisSize < size) {	size = thisSize;	}
				}
			}else { //make sure size is not too large
				vector<string> newGroups;
				for (int i = 0; i < Groups.size(); i++) {
					int thisSize = groupMap->getNumSeqs(Groups[i]);
					
					if (thisSize >= size) {	newGroups.push_back(Groups[i]);	}
					else {  m->mothurOut("You have selected a size that is larger than " + Groups[i] + " number of sequences, removing " + Groups[i] + "."); m->mothurOutEndLine(); }
				}
				Groups = newGroups;
			}
			
			m->mothurOut("Sampling " + toString(size) + " from each group."); m->mothurOutEndLine();	
		}else{
			if (pickedGroups) {
				int total = 0;
				for(int i = 0; i < Groups.size(); i++) {
					total += groupMap->getNumSeqs(Groups[i]);
				}
				
				if (size == 0) { //user has not set size, set size = 10% samples size
					size = int (total * 0.10);
				}
				
				if (total < size) { 
					m->mothurOut("Your size is too large for the number of groups you selected. Adjusting to " + toString(int (total * 0.10)) + "."); m->mothurOutEndLine();
					size = int (total * 0.10);
				}
				
				m->mothurOut("Sampling " + toString(size) + " from " + toString(total) + "."); m->mothurOutEndLine();
			}else{
				
				if (size == 0) { //user has not set size, set size = 10% samples size
					size = int (list->getNumSeqs() * 0.10);
				}
				
				int thisSize = list->getNumSeqs();
				if (size > thisSize) { m->mothurOut("Your list file only contains " + toString(thisSize) + " sequences. Setting size to " + toString(thisSize) + "."); m->mothurOutEndLine();
					size = thisSize;
				}
				
				m->mothurOut("Sampling " + toString(size) + " from " + toString(list->getNumSeqs()) + "."); m->mothurOutEndLine();
			}
		}
		
		
		//fill names
		for (int i = 0; i < list->getNumBins(); i++) {
			string binnames = list->get(i);
			
			//parse names
			string individual = "";
			int length = binnames.length();
			for(int j=0;j<length;j++){
				if(binnames[j] == ','){
					
					if (groupfile != "") { //if there is a groupfile given fill in group info
						string group = groupMap->getGroup(individual);
						if (group == "not found") { m->mothurOut("[ERROR]: " + individual + " is not in your groupfile. please correct."); m->mothurOutEndLine(); group = "NOTFOUND"; }
						
						if (pickedGroups) { //if hte user picked groups, we only want to keep the names of sequences from those groups
							if (m->inUsersGroups(group, Groups)) {
								names.push_back(individual);
							}
						}else{
							names.push_back(individual);
						}
					}else{ //save everyone, group
						names.push_back(individual);
					}
					individual = "";				
				}
				else{
					individual += binnames[j];
				}
			}
			//save last name
			if (groupfile != "") { //if there is a groupfile given fill in group info
				string group = groupMap->getGroup(individual);
				if (group == "not found") { m->mothurOut("[ERROR]: " + individual + " is not in your groupfile. please correct."); m->mothurOutEndLine(); group = "NOTFOUND"; }
				
				if (pickedGroups) { //if hte user picked groups, we only want to keep the names of sequences from those groups
					if (m->inUsersGroups(group, Groups)) {
						names.push_back(individual);
					}
				}else{
					names.push_back(individual);
				}
			}else{ //save everyone, group
				names.push_back(individual);
			}
		}
		
		random_shuffle(names.begin(), names.end());
			
		//randomly select a subset of those names to include in the subsample
		set<string> subset; //dont want repeat sequence names added
		if (persample) {
			for (int i = 0; i < Groups.size(); i++) {
				
				for (int j = 0; j < size; j++) {
					
					if (m->control_pressed) { break; }
					
					//get random sequence to add, making sure we have not already added it
					bool done = false;
					int myrand;
					while (!done) {
						myrand = int((float)(names.size()) * (float)(rand()) / ((float)RAND_MAX+1.0));
						
						if (subset.count(names[myrand]) == 0) { //you are not already added
							if (groupMap->getGroup(names[myrand]) == Groups[i])  { subset.insert(names[myrand]); break;	}
						}
					}
				}
			}
		}else{
			for (int j = 0; j < size; j++) {
				
				if (m->control_pressed) { break; }
				
				//get random sequence to add, making sure we have not already added it
				bool done = false;
				int myrand;
				while (!done) {
					myrand = int((float)(names.size()) * (float)(rand()) / ((float)RAND_MAX+1.0));
					
					if (subset.count(names[myrand]) == 0)  { subset.insert(names[myrand]); break;	}
				}
			}	
		}
		
		if (groupfile != "") { 
			//write out new groupfile
			for (set<string>::iterator it = subset.begin(); it != subset.end(); it++) {
				string group = groupMap->getGroup(*it);
				if (group == "not found") { group = "NOTFOUND"; }
				
				outGroup << *it << '\t' << group << endl;
			}
			outGroup.close(); delete groupMap; 
		}
		
						
		//as long as you are not at the end of the file or done wih the lines you want
		while((list != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			
			if (m->control_pressed) {  delete list; delete input; out.close();  return 0;  }
			
			if(allLines == 1 || labels.count(list->getLabel()) == 1){			
				
				m->mothurOut(list->getLabel()); m->mothurOutEndLine();
				
				processList(list, out, subset);
				
				processedLabels.insert(list->getLabel());
				userLabels.erase(list->getLabel());
			}
			
			if ((m->anyLabelsToProcess(list->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
				string saveLabel = list->getLabel();
				
				delete list; 
				
				list = input->getListVector(lastLabel);
				m->mothurOut(list->getLabel()); m->mothurOutEndLine();
				
				processList(list, out, subset);
				
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
		
		
		if (m->control_pressed) {  if (list != NULL) { delete list; } delete input; out.close(); return 0;  }
		
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
			
			processList(list, out, subset);
			
			delete list; list = NULL;
		}
		
		out.close();  
		if (list != NULL) { delete list; }
		delete input;
						
		return 0;
 
	}
	catch(exception& e) {
		m->errorOut(e, "SubSampleCommand", "getSubSampleList");
		exit(1);
	}
}
//**********************************************************************************************************************
int SubSampleCommand::processList(ListVector*& list, ofstream& out, set<string>& subset) {
	try {
				
		int numBins = list->getNumBins();

		ListVector* temp = new ListVector();
		temp->setLabel(list->getLabel());
		
		for (int i = 0; i < numBins; i++) {
			
			if (m->control_pressed) { break; }
			
			string binnames = list->get(i);
			
			//parse names
			string individual = "";
			string newNames = "";
			int length = binnames.length();
			for(int j=0;j<length;j++){
				if(binnames[j] == ','){
					if (subset.count(individual) != 0) {  newNames += individual + ",";  }
					individual = "";				
				}else{
					individual += binnames[j];
				}
			}
			if (subset.count(individual) != 0) {  newNames += individual + ",";  }
			
			
			//if there are names in this bin add to new list
			if (newNames != "") { 
				newNames = newNames.substr(0, newNames.length()-1); //rip off extra comma
				temp->push_back(newNames);
			}
		}
		
		delete list;
		list = temp;
		
		if (m->control_pressed) { return 0; }
		
		list->print(out);
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "SubSampleCommand", "processList");
		exit(1);
	}
}
//**********************************************************************************************************************
int SubSampleCommand::getSubSampleRabund() {
	try {
		InputData* input = new InputData(rabundfile, "rabund");
		RAbundVector* rabund = input->getRAbundVector();
		string lastLabel = rabund->getLabel();
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;
		
		if (size == 0) { //user has not set size, set size = 10%
			size = int((rabund->getNumSeqs()) * 0.10);
		}else if (size > rabund->getNumSeqs()) { m->mothurOut("The size you selected is too large, skipping rabund file."); m->mothurOutEndLine(); delete input; delete rabund; return 0; }
		
		m->mothurOut("Sampling " + toString(size) + " from " + toString(rabund->getNumSeqs()) + "."); m->mothurOutEndLine();
		
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += m->hasPath(rabundfile);  }
		string outputFileName = thisOutputDir + m->getRootName(m->getSimpleName(rabundfile)) + "subsample" + m->getExtension(rabundfile);
		
		ofstream out;
		m->openOutputFile(outputFileName, out);
		outputTypes["rabund"].push_back(outputFileName);  outputNames.push_back(outputFileName);
		
		//as long as you are not at the end of the file or done wih the lines you want
		while((rabund != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			if (m->control_pressed) {  delete input; delete rabund; out.close(); return 0;  }
			
			if(allLines == 1 || labels.count(rabund->getLabel()) == 1){			
				
				m->mothurOut(rabund->getLabel()); m->mothurOutEndLine();
				
				processRabund(rabund, out);
				
				processedLabels.insert(rabund->getLabel());
				userLabels.erase(rabund->getLabel());
			}
			
			if ((m->anyLabelsToProcess(rabund->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
				string saveLabel = rabund->getLabel();
				
				delete rabund; 
				
				rabund = input->getRAbundVector(lastLabel);
				m->mothurOut(rabund->getLabel()); m->mothurOutEndLine();
				
				processRabund(rabund, out);
				
				processedLabels.insert(rabund->getLabel());
				userLabels.erase(rabund->getLabel());
				
				//restore real lastlabel to save below
				rabund->setLabel(saveLabel);
			}
			
			lastLabel = rabund->getLabel();
			
			//prevent memory leak
			delete rabund; rabund = NULL;
			
			//get next line to process
			rabund = input->getRAbundVector();				
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
			if (rabund != NULL) { delete rabund; }
			
			rabund = input->getRAbundVector(lastLabel);
			
			m->mothurOut(rabund->getLabel()); m->mothurOutEndLine();
			
			processRabund(rabund, out);
			
			delete rabund;
		}
		
		delete input;
		out.close();  
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "SubSampleCommand", "getSubSampleRabund");
		exit(1);
	}
}
//**********************************************************************************************************************
int SubSampleCommand::processRabund(RAbundVector*& rabund, ofstream& out) {
	try {
		
		int numBins = rabund->getNumBins();
		int thisSize = rabund->getNumSeqs();
			
		if (thisSize != size) {
				
			OrderVector* order = new OrderVector();
			for(int p=0;p<numBins;p++){
				for(int j=0;j<rabund->get(p);j++){
					order->push_back(p);
				}
			}
			random_shuffle(order->begin(), order->end());
			
			RAbundVector* temp = new RAbundVector(numBins);
			temp->setLabel(rabund->getLabel());
			
			delete rabund;
			rabund = temp;
			
			for (int j = 0; j < size; j++) {
				
				if (m->control_pressed) { delete order; return 0; }
				
				//get random number to sample from order between 0 and thisSize-1.
				int myrand = int((float)(thisSize) * (float)(rand()) / ((float)RAND_MAX+1.0));
				
				int bin = order->get(myrand);
				
				int abund = rabund->get(bin);
				rabund->set(bin, (abund+1));
			}
			
			delete order;
		}
		
		if (m->control_pressed) { return 0; }
		
		rabund->print(out);
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "SubSampleCommand", "processRabund");
		exit(1);
	}
}	
//**********************************************************************************************************************
int SubSampleCommand::getSubSampleSabund() {
	try {
				
		InputData* input = new InputData(sabundfile, "sabund");
		SAbundVector* sabund = input->getSAbundVector();
		string lastLabel = sabund->getLabel();
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;
		
		if (size == 0) { //user has not set size, set size = 10%
			size = int((sabund->getNumSeqs()) * 0.10);
		}else if (size > sabund->getNumSeqs()) { m->mothurOut("The size you selected is too large, skipping sabund file."); m->mothurOutEndLine(); delete input; delete sabund; return 0; }
		
		
		m->mothurOut("Sampling " + toString(size) + " from " + toString(sabund->getNumSeqs()) + "."); m->mothurOutEndLine();
		
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += m->hasPath(sabundfile);  }
		string outputFileName = thisOutputDir + m->getRootName(m->getSimpleName(sabundfile)) + "subsample" + m->getExtension(sabundfile);
		
		ofstream out;
		m->openOutputFile(outputFileName, out);
		outputTypes["sabund"].push_back(outputFileName);  outputNames.push_back(outputFileName);
		
		
		//as long as you are not at the end of the file or done wih the lines you want
		while((sabund != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			if (m->control_pressed) {  delete input; delete sabund; out.close(); return 0;  }
			
			if(allLines == 1 || labels.count(sabund->getLabel()) == 1){			
				
				m->mothurOut(sabund->getLabel()); m->mothurOutEndLine();
				
				processSabund(sabund, out);
				
				processedLabels.insert(sabund->getLabel());
				userLabels.erase(sabund->getLabel());
			}
			
			if ((m->anyLabelsToProcess(sabund->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
				string saveLabel = sabund->getLabel();
				
				delete sabund; 
				
				sabund = input->getSAbundVector(lastLabel);
				m->mothurOut(sabund->getLabel()); m->mothurOutEndLine();
				
				processSabund(sabund, out);
				
				processedLabels.insert(sabund->getLabel());
				userLabels.erase(sabund->getLabel());
				
				//restore real lastlabel to save below
				sabund->setLabel(saveLabel);
			}
			
			lastLabel = sabund->getLabel();
			
			//prevent memory leak
			delete sabund; sabund = NULL;
			
			//get next line to process
			sabund = input->getSAbundVector();				
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
			if (sabund != NULL) { delete sabund; }
			
			sabund = input->getSAbundVector(lastLabel);
			
			m->mothurOut(sabund->getLabel()); m->mothurOutEndLine();
			
			processSabund(sabund, out);
			
			delete sabund;
		}
		
		delete input;
		out.close();  
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "SubSampleCommand", "getSubSampleSabund");
		exit(1);
	}
}
//**********************************************************************************************************************
int SubSampleCommand::processSabund(SAbundVector*& sabund, ofstream& out) {
	try {
		
		RAbundVector* rabund = new RAbundVector();
		*rabund = sabund->getRAbundVector();
		
		int numBins = rabund->getNumBins();
		int thisSize = rabund->getNumSeqs();
	
		if (thisSize != size) {
			
			OrderVector* order = new OrderVector();
			for(int p=0;p<numBins;p++){
				for(int j=0;j<rabund->get(p);j++){
					order->push_back(p);
				}
			}
			random_shuffle(order->begin(), order->end());
			
			RAbundVector* temp = new RAbundVector(numBins);
			temp->setLabel(rabund->getLabel());
			
			delete rabund;
			rabund = temp;
			
			for (int j = 0; j < size; j++) {
	
				if (m->control_pressed) { delete order; return 0; }
				
				//get random number to sample from order between 0 and thisSize-1.
				int myrand = int((float)(thisSize) * (float)(rand()) / ((float)RAND_MAX+1.0));
				
				int bin = order->get(myrand);
				
				int abund = rabund->get(bin);
				rabund->set(bin, (abund+1));
			}
			
			delete order;
		}
		
		if (m->control_pressed) { return 0; }

		delete sabund;
		sabund = new SAbundVector();
		*sabund = rabund->getSAbundVector();
		delete rabund;
	
		sabund->print(out);
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "SubSampleCommand", "processSabund");
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



