/*
 *  getsharedotucommand.cpp
 *  Mothur
 *
 *  Created by westcott on 9/22/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */

#include "getsharedotucommand.h"


//**********************************************************************************************************************
vector<string> GetSharedOTUCommand::getValidParameters(){	
	try {
		string Array[] =  {"label","unique","shared","fasta","list","group","output","outputdir","inputdir"};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "GetSharedOTUCommand", "getValidParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
GetSharedOTUCommand::GetSharedOTUCommand(){	
	try {
		//initialize outputTypes
		vector<string> tempOutNames;
		outputTypes["fasta"] = tempOutNames;
		outputTypes["accnos"] = tempOutNames;
		outputTypes["sharedseqs"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "GetSharedOTUCommand", "GetSharedOTUCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> GetSharedOTUCommand::getRequiredParameters(){	
	try {
		string Array[] =  {"list","group"};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "GetSharedOTUCommand", "getRequiredParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> GetSharedOTUCommand::getRequiredFiles(){	
	try {
		vector<string> myArray;
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "GetSharedOTUCommand", "getRequiredFiles");
		exit(1);
	}
}
//**********************************************************************************************************************
GetSharedOTUCommand::GetSharedOTUCommand(string option)  {
	try {
	
		globaldata = GlobalData::getInstance();
		abort = false;
		unique = true;
		allLines = 1;
		labels.clear();
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"label","unique","shared","fasta","list","group","output","outputdir","inputdir"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			map<string,string>::iterator it;
			
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["fasta"] = tempOutNames;
			outputTypes["accnos"] = tempOutNames;
			outputTypes["sharedseqs"] = tempOutNames;
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";		}
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("fasta");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["fasta"] = inputDir + it->second;		}
				}
				
				it = parameters.find("list");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["list"] = inputDir + it->second;		}
				}
				
				it = parameters.find("group");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["group"] = inputDir + it->second;		}
				}
			}

			
			//check for required parameters
			listfile = validParameter.validFile(parameters, "list", true);
			if (listfile == "not open") { abort = true; }
			else if (listfile == "not found") { listfile = ""; }	
			else {  globaldata->setListFile(listfile);  globaldata->setFormat("list"); 	}
			
			groupfile = validParameter.validFile(parameters, "group", true);
			if (groupfile == "not open") { abort = true; }	
			else if (groupfile == "not found") { groupfile = ""; }
						
			if ((listfile == "") || (groupfile == "")) { m->mothurOut("The list and group parameters are required."); m->mothurOutEndLine(); abort = true; }
			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  m->splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
			
			output = validParameter.validFile(parameters, "output", false);			
			if (output == "not found") { output = ""; }
			
			groups = validParameter.validFile(parameters, "unique", false);			
			if (groups == "not found") { groups = ""; }
			else { 
				userGroups = "unique." + groups;
				m->splitAtDash(groups, Groups);
				globaldata->Groups = Groups;
				
			}
			
			groups = validParameter.validFile(parameters, "shared", false);			
			if (groups == "not found") { groups = "";  }
			else { 
				userGroups = groups;
				m->splitAtDash(groups, Groups);
				globaldata->Groups = Groups;
				unique = false;
			}
			
			fastafile = validParameter.validFile(parameters, "fasta", true);
			if (fastafile == "not open") { abort = true; }
			else if (fastafile == "not found") {  fastafile = "";  }	
				
		}

	}
	catch(exception& e) {
		m->errorOut(e, "GetSharedOTUCommand", "GetSharedOTUCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

void GetSharedOTUCommand::help(){
	try {
		m->mothurOut("The get.sharedseqs command parameters are list, group, label, unique, shared, output and fasta.  The list and group parameters are required.\n");
		m->mothurOut("The label parameter allows you to select what distance levels you would like output files for, and are separated by dashes.\n");
		m->mothurOut("The unique and shared parameters allow you to select groups you would like to know the shared info for, and are separated by dashes.\n");
		m->mothurOut("If you enter your groups under the unique parameter mothur will return the otus that contain ONLY sequences from those groups.\n");
		m->mothurOut("If you enter your groups under the shared parameter mothur will return the otus that contain sequences from those groups and may also contain sequences from other groups.\n");
		m->mothurOut("If you do not enter any groups then the get.sharedseqs command will return sequences that are unique to all groups in your group file.\n");
		m->mothurOut("The fasta parameter allows you to input a fasta file and outputs a fasta file for each distance level containing only the sequences that are in OTUs shared by the groups specified.\n");
		m->mothurOut("The output parameter allows you to output the list of names without the group and bin number added. \n");
		m->mothurOut("With this option you can use the names file as an input in get.seqs and remove.seqs commands. To do this enter output=accnos. \n");
		m->mothurOut("The get.sharedseqs command outputs a .names file for each distance level containing a list of sequences in the OTUs shared by the groups specified.\n");
		m->mothurOut("The get.sharedseqs command should be in the following format: get.sabund(label=yourLabels, groups=yourGroups, fasta=yourFastafile, output=yourOutput).\n");
		m->mothurOut("Example get.sharedseqs(list=amazon.fn.list, label=unique-0.01, group=forest-pasture, fasta=amazon.fasta, output=accnos).\n");
		m->mothurOut("The output to the screen is the distance and the number of otus at that distance for the groups you specified.\n");
		m->mothurOut("The default value for label is all labels in your inputfile. The default for groups is all groups in your file.\n");
		m->mothurOut("Note: No spaces between parameter labels (i.e. label), '=' and parameters (i.e.yourLabel).\n\n");
	}
	catch(exception& e) {
		m->errorOut(e, "GetSharedOTUCommand", "help");
		exit(1);
	}
}

//**********************************************************************************************************************

GetSharedOTUCommand::~GetSharedOTUCommand(){}

//**********************************************************************************************************************

int GetSharedOTUCommand::execute(){
	try {
		
		if (abort == true) { return 0; }
		
		groupMap = new GroupMap(groupfile);
		int error = groupMap->readMap();
		if (error == 1) { delete groupMap; return 0; }
		
		if (m->control_pressed) { delete groupMap; return 0; }
		
		globaldata->gGroupmap = groupMap;
		
		if (Groups.size() == 0) {
			Groups = groupMap->namesOfGroups;
			
			//make string for outputfile name
			userGroups = "unique.";
			for(int i = 0; i < Groups.size(); i++) {  userGroups += Groups[i] + "-";  }
			userGroups = userGroups.substr(0, userGroups.length()-1);
		}
	
		//put groups in map to find easier
		for(int i = 0; i < Groups.size(); i++) {
			groupFinder[Groups[i]] = Groups[i];
		}
		
		if (fastafile != "") {
			ifstream inFasta;
			m->openInputFile(fastafile, inFasta);
			
			while(!inFasta.eof()) {
				if (m->control_pressed) { outputTypes.clear(); inFasta.close(); delete groupMap; return 0; }
				
				Sequence seq(inFasta); m->gobble(inFasta);
				if (seq.getName() != "") {  seqs.push_back(seq);   }
			}
			inFasta.close();
		}
		
		ListVector* lastlist = NULL;
		string lastLabel = "";
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;
		
		ifstream in;
		m->openInputFile(listfile, in);
		
		//as long as you are not at the end of the file or done wih the lines you want
		while((!in.eof()) && ((allLines == 1) || (userLabels.size() != 0))) {
			
			if (m->control_pressed) { 
				if (lastlist != NULL) {		delete lastlist;	}
				for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); }  outputTypes.clear();
				delete groupMap; return 0;
			}
			
			list = new ListVector(in);
			
			if(allLines == 1 || labels.count(list->getLabel()) == 1){			
				m->mothurOut(list->getLabel()); 
				process(list);
				
				processedLabels.insert(list->getLabel());
				userLabels.erase(list->getLabel());
			}
			
			if ((m->anyLabelsToProcess(list->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
					string saveLabel = list->getLabel();
					
					m->mothurOut(lastlist->getLabel()); 
					process(lastlist);
					
					processedLabels.insert(lastlist->getLabel());
					userLabels.erase(lastlist->getLabel());
					
					//restore real lastlabel to save below
					list->setLabel(saveLabel);
			}

			lastLabel = list->getLabel();
			
			if (lastlist != NULL) {		delete lastlist;	}
			lastlist = list;			
		}
		
		in.close();
		
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
				m->mothurOut(lastlist->getLabel()); 
				process(lastlist);
					
				processedLabels.insert(lastlist->getLabel());
				userLabels.erase(lastlist->getLabel());
		}
		

		//reset groups parameter
		globaldata->Groups.clear();  
		
		if (lastlist != NULL) {		delete lastlist;	}
		
		if (m->control_pressed) { outputTypes.clear(); for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); }  delete groupMap; return 0; } 
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();


		return 0;
	}

	catch(exception& e) {
		m->errorOut(e, "GetSharedOTUCommand", "execute");
		exit(1);
	}
}
/***********************************************************/
int GetSharedOTUCommand::process(ListVector* shared) {
	try {
		
		map<string, string> fastaMap;
		
		ofstream outNames;
		string outputFileNames;
		
		if (outputDir == "") { outputDir += m->hasPath(listfile); }
		if (output != "accnos") {
			outputFileNames = outputDir + m->getRootName(m->getSimpleName(listfile)) + shared->getLabel() + userGroups + ".shared.seqs";
		}else {
			outputFileNames = outputDir + m->getRootName(m->getSimpleName(listfile)) + shared->getLabel() + userGroups + ".accnos";
		}
		m->openOutputFile(outputFileNames, outNames);
		
		bool wroteSomething = false;
		int num = 0;
				
		//go through each bin, find out if shared
		for (int i = 0; i < shared->getNumBins(); i++) {
			if (m->control_pressed) { outNames.close(); remove(outputFileNames.c_str()); return 0; }
			
			bool uniqueOTU = true;
			
			map<string, int> atLeastOne;
			for (int f = 0; f < Groups.size(); f++) {
				atLeastOne[Groups[f]] = 0;
			}
			
			vector<string> namesOfSeqsInThisBin;
			
			string names = shared->get(i);  
			while ((names.find_first_of(',') != -1)) { 
				string name = names.substr(0,names.find_first_of(','));
				names = names.substr(names.find_first_of(',')+1, names.length());
				
				//find group
				string seqGroup = groupMap->getGroup(name);
				if (output != "accnos") {
					namesOfSeqsInThisBin.push_back((name + "|" + seqGroup + "|" + toString(i+1)));
				}else {  namesOfSeqsInThisBin.push_back(name);	}
				
				if (seqGroup == "not found") { m->mothurOut(name + " is not in your groupfile. Please correct."); m->mothurOutEndLine(); exit(1);  }
				
				//is this seq in one of hte groups we care about
				it = groupFinder.find(seqGroup);
				if (it == groupFinder.end()) {  uniqueOTU = false;  } //you have a sequence from a group you don't want
				else {  atLeastOne[seqGroup]++;  }
			}
			
			//get last name
			string seqGroup = groupMap->getGroup(names);
			if (output != "accnos") {
				namesOfSeqsInThisBin.push_back((names + "|" + seqGroup + "|" + toString(i+1)));
			}else {  namesOfSeqsInThisBin.push_back(names);	}
			
			if (seqGroup == "not found") { m->mothurOut(names + " is not in your groupfile. Please correct."); m->mothurOutEndLine(); exit(1);  }
			
			//is this seq in one of hte groups we care about
			it = groupFinder.find(seqGroup);
			if (it == groupFinder.end()) {  uniqueOTU = false;  } //you have a sequence from a group you don't want
			else {  atLeastOne[seqGroup]++;  }
			
			
			//make sure you have at least one seq from each group you want
			bool sharedByAll = true;
			map<string, int>::iterator it2;
			for (it2 = atLeastOne.begin(); it2 != atLeastOne.end(); it2++) {
				if (it2->second == 0) {  sharedByAll = false;	}
			}
			
			//if the user wants unique bins and this is unique then print
			//or this the user wants shared bins and this bin is shared then print
			if ((unique && uniqueOTU && sharedByAll) || (!unique && sharedByAll)) {
				
				wroteSomething = true;
				num++;
				
				//output list of names 
				for (int j = 0; j < namesOfSeqsInThisBin.size(); j++) {
					outNames << namesOfSeqsInThisBin[j] << endl;
					
					if (fastafile != "") { 
						if (output != "accnos") {
							string seqName = namesOfSeqsInThisBin[j].substr(0,namesOfSeqsInThisBin[j].find_last_of('|'));
							seqName = seqName.substr(0,seqName.find_last_of('|'));
							fastaMap[seqName] = namesOfSeqsInThisBin[j];  //fastaMap needs to contain just the seq name for output later
						}else {
							fastaMap[namesOfSeqsInThisBin[j]] = namesOfSeqsInThisBin[j];
						}
					}
				}
			}
		}
		
		outNames.close();
		
		if (!wroteSomething) {
			remove(outputFileNames.c_str());
			string outputString = "\t" + toString(num) + " - No otus shared by groups";
			
			string groupString = "";
			for (int h = 0; h < Groups.size(); h++) {
				groupString += "  " + Groups[h];
			}
			
			outputString += groupString + ".";
			m->mothurOut(outputString); m->mothurOutEndLine();
		}else { 
			m->mothurOut("\t" + toString(num)); m->mothurOutEndLine(); 
			outputNames.push_back(outputFileNames);
			if (output != "accnos") { outputTypes["sharedseqs"].push_back(outputFileNames); }
			else { outputTypes["accnos"].push_back(outputFileNames); }
		}
		
		//if fasta file provided output new fasta file
		if ((fastafile != "") && wroteSomething) {
			if (outputDir == "") { outputDir += m->hasPath(fastafile); }
			string outputFileFasta = outputDir + m->getRootName(m->getSimpleName(fastafile)) + shared->getLabel() + userGroups + ".shared.fasta";
			ofstream outFasta;
			m->openOutputFile(outputFileFasta, outFasta);
			outputNames.push_back(outputFileFasta); outputTypes["fasta"].push_back(outputFileFasta);
			
			for (int k = 0; k < seqs.size(); k++) {
				if (m->control_pressed) { outFasta.close(); return 0; }
			
				//if this is a sequence we want, output it
				it = fastaMap.find(seqs[k].getName());
				if (it != fastaMap.end()) {
				
					if (output != "accnos") {
						outFasta << ">" << it->second << endl;
					}else {
						outFasta << ">" << it->first << endl;
					}
					
					outFasta << seqs[k].getAligned() << endl;
				}
			}
			
			outFasta.close();
		}
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "GetSharedOTUCommand", "process");
		exit(1);
	}
}

//**********************************************************************************************************************
