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

GetSharedOTUCommand::GetSharedOTUCommand(string option){
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
					path = hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["fasta"] = inputDir + it->second;		}
				}
				
				it = parameters.find("list");
				//user has given a template file
				if(it != parameters.end()){ 
					path = hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["list"] = inputDir + it->second;		}
				}
				
				it = parameters.find("group");
				//user has given a template file
				if(it != parameters.end()){ 
					path = hasPath(it->second);
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
						
			if ((listfile == "") || (groupfile == "")) { mothurOut("The list and group parameters are required."); mothurOutEndLine(); abort = true; }
			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
			
			output = validParameter.validFile(parameters, "output", false);			
			if (output == "not found") { output = ""; }
			
			groups = validParameter.validFile(parameters, "unique", false);			
			if (groups == "not found") { groups = ""; }
			else { 
				userGroups = "unique." + groups;
				splitAtDash(groups, Groups);
				globaldata->Groups = Groups;
				
			}
			
			groups = validParameter.validFile(parameters, "shared", false);			
			if (groups == "not found") { groups = "";  }
			else { 
				userGroups = groups;
				splitAtDash(groups, Groups);
				globaldata->Groups = Groups;
				unique = false;
			}
			
			fastafile = validParameter.validFile(parameters, "fasta", true);
			if (fastafile == "not open") { abort = true; }
			else if (fastafile == "not found") {  fastafile = "";  }	
				
		}

	}
	catch(exception& e) {
		errorOut(e, "GetSharedOTUCommand", "GetSharedOTUCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

void GetSharedOTUCommand::help(){
	try {
		mothurOut("The get.sharedseqs command parameters are list, group, label, unique, shared, output and fasta.  The list and group parameters are required.\n");
		mothurOut("The label parameter allows you to select what distance levels you would like output files for, and are separated by dashes.\n");
		mothurOut("The unique and shared parameters allow you to select groups you would like to know the shared info for, and are separated by dashes.\n");
		mothurOut("If you enter your groups under the unique parameter mothur will return the otus that contain ONLY sequences from those groups.\n");
		mothurOut("If you enter your groups under the shared parameter mothur will return the otus that contain sequences from those groups and may also contain sequences from other groups.\n");
		mothurOut("If you do not enter any groups then the get.sharedseqs command will return sequences that are unique to all groups in your group file.\n");
		mothurOut("The fasta parameter allows you to input a fasta file and outputs a fasta file for each distance level containing only the sequences that are in OTUs shared by the groups specified.\n");
		mothurOut("The output parameter allows you to output the list of names without the group and bin number added. \n");
		mothurOut("With this option you can use the names file as an input in get.seqs and remove.seqs commands. To do this enter output=accnos. \n");
		mothurOut("The get.sharedseqs command outputs a .names file for each distance level containing a list of sequences in the OTUs shared by the groups specified.\n");
		mothurOut("The get.sharedseqs command should be in the following format: get.sabund(label=yourLabels, groups=yourGroups, fasta=yourFastafile, output=yourOutput).\n");
		mothurOut("Example get.sharedseqs(list=amazon.fn.list, label=unique-0.01, group=forest-pasture, fasta=amazon.fasta, output=accnos).\n");
		mothurOut("The output to the screen is the distance and the number of otus at that distance for the groups you specified.\n");
		mothurOut("The default value for label is all labels in your inputfile. The default for groups is all groups in your file.\n");
		mothurOut("Note: No spaces between parameter labels (i.e. label), '=' and parameters (i.e.yourLabel).\n\n");
	}
	catch(exception& e) {
		errorOut(e, "GetSharedOTUCommand", "help");
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
			openInputFile(fastafile, inFasta);
			
			while(!inFasta.eof()) {
				Sequence seq(inFasta); gobble(inFasta);
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
		openInputFile(listfile, in);
		
		//as long as you are not at the end of the file or done wih the lines you want
		while((!in.eof()) && ((allLines == 1) || (userLabels.size() != 0))) {
			
			list = new ListVector(in);
			
			if(allLines == 1 || labels.count(list->getLabel()) == 1){			
				mothurOut(list->getLabel()); 
				process(list);
				
				processedLabels.insert(list->getLabel());
				userLabels.erase(list->getLabel());
			}
			
			if ((anyLabelsToProcess(list->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
					string saveLabel = list->getLabel();
					
					mothurOut(lastlist->getLabel()); 
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
				mothurOut(lastlist->getLabel()); 
				process(lastlist);
					
				processedLabels.insert(lastlist->getLabel());
				userLabels.erase(lastlist->getLabel());
		}
		

		//reset groups parameter
		globaldata->Groups.clear();  
		
		if (lastlist != NULL) {		delete lastlist;	}
		return 0;
	}

	catch(exception& e) {
		errorOut(e, "GetSharedOTUCommand", "execute");
		exit(1);
	}
}
/***********************************************************/
void GetSharedOTUCommand::process(ListVector* shared) {
	try {
		
		map<string, string> fastaMap;
		
		ofstream outNames;
		string outputFileNames;
		
		if (outputDir == "") { outputDir += hasPath(listfile); }
		if (output != "accnos") {
			outputFileNames = outputDir + getRootName(getSimpleName(listfile)) + shared->getLabel() + userGroups + ".shared.seqs";
		}else {
			outputFileNames = outputDir + getRootName(getSimpleName(listfile)) + shared->getLabel() + userGroups + ".accnos";
		}
		openOutputFile(outputFileNames, outNames);
		
		bool wroteSomething = false;
		int num = 0;
				
		//go through each bin, find out if shared
		for (int i = 0; i < shared->getNumBins(); i++) {
			
			bool uniqueOTU = true;
			
			map<string, int> atLeastOne;
			for (int m = 0; m < Groups.size(); m++) {
				atLeastOne[Groups[m]] = 0;
			}
			
			vector<string> namesOfSeqsInThisBin;
			
			string names = shared->get(i);  
			while ((names.find_first_of(',') != -1)) { 
				string name = names.substr(0,names.find_first_of(','));
				names = names.substr(names.find_first_of(',')+1, names.length());
				
				//find group
				string seqGroup = groupMap->getGroup(name);
				if (output != "accnos") {
					namesOfSeqsInThisBin.push_back((name + "\t" + seqGroup + "\t" + toString(i+1)));
				}else {  namesOfSeqsInThisBin.push_back(name);	}
				
				if (seqGroup == "not found") { mothurOut(name + " is not in your groupfile. Please correct."); mothurOutEndLine(); exit(1);  }
				
				//is this seq in one of hte groups we care about
				it = groupFinder.find(seqGroup);
				if (it == groupFinder.end()) {  uniqueOTU = false;  } //you have a sequence from a group you don't want
				else {  atLeastOne[seqGroup]++;  }
			}
			
			//get last name
			string seqGroup = groupMap->getGroup(names);
			if (output != "accnos") {
				namesOfSeqsInThisBin.push_back((names + "\t" + seqGroup + "\t" + toString(i+1)));
			}else {  namesOfSeqsInThisBin.push_back(names);	}
			
			if (seqGroup == "not found") { mothurOut(names + " is not in your groupfile. Please correct."); mothurOutEndLine(); exit(1);  }
			
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
							string seqName = namesOfSeqsInThisBin[j].substr(0,namesOfSeqsInThisBin[j].find_last_of('\t'));
							seqName = seqName.substr(0,seqName.find_last_of('\t'));
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
			mothurOut(outputString); mothurOutEndLine();
		}else { mothurOut("\t" + toString(num)); mothurOutEndLine(); }
		
		//if fasta file provided output new fasta file
		if ((fastafile != "") && wroteSomething) {
			if (outputDir == "") { outputDir += hasPath(fastafile); }
			string outputFileFasta = outputDir + getRootName(getSimpleName(fastafile)) + shared->getLabel() + userGroups + ".shared.fasta";
			ofstream outFasta;
			openOutputFile(outputFileFasta, outFasta);
			
			for (int k = 0; k < seqs.size(); k++) {
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

		
	}
	catch(exception& e) {
		errorOut(e, "GetSharedOTUCommand", "process");
		exit(1);
	}
}

//**********************************************************************************************************************
