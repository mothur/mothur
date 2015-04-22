/*
 *  getsharedotucommand.cpp
 *  Mothur
 *
 *  Created by westcott on 9/22/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */

#include "getsharedotucommand.h"
#include "sharedutilities.h"

//**********************************************************************************************************************
vector<string> GetSharedOTUCommand::setParameters(){	
	try {
		CommandParameter pfasta("fasta", "InputTypes", "", "", "sharedFasta", "none", "none","fasta",false,false); parameters.push_back(pfasta);
		CommandParameter pgroup("group", "InputTypes", "", "", "none", "none", "groupList","",false,false,true); parameters.push_back(pgroup);
		CommandParameter plist("list", "InputTypes", "", "", "sharedList", "sharedList", "groupList","sharedseq",false,false,true); parameters.push_back(plist);
        CommandParameter pshared("shared", "InputTypes", "", "", "sharedList-sharedFasta", "sharedList", "none","sharedseq",false,false,true); parameters.push_back(pshared);
		CommandParameter poutput("output", "Multiple", "accnos-default", "default", "", "", "","",false,false); parameters.push_back(poutput);
		CommandParameter plabel("label", "String", "", "", "", "", "","",false,false); parameters.push_back(plabel);
		CommandParameter puniquegroups("uniquegroups", "String", "", "", "", "", "","",false,false,true); parameters.push_back(puniquegroups);
		CommandParameter psharedgroups("sharedgroups", "String", "", "", "", "", "","",false,false,true); parameters.push_back(psharedgroups);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);

		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "GetSharedOTUCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string GetSharedOTUCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The get.sharedseqs command parameters are list, group, shared, label, uniquegroups, sharedgroups, output and fasta.  The list and group or shared parameters are required, unless you have valid current files.\n";
		helpString += "The label parameter allows you to select what distance levels you would like output files for, and are separated by dashes.\n";
		helpString += "The uniquegroups and sharedgroups parameters allow you to select groups you would like to know the shared info for, and are separated by dashes.\n";
		helpString += "If you enter your groups under the uniquegroups parameter mothur will return the otus that contain ONLY sequences from those groups.\n";
		helpString += "If you enter your groups under the sharedgroups parameter mothur will return the otus that contain sequences from those groups and may also contain sequences from other groups.\n";
		helpString += "If you do not enter any groups then the get.sharedseqs command will return sequences that are unique to all groups in your group or shared file.\n";
		helpString += "The fasta parameter allows you to input a fasta file and outputs a fasta file for each distance level containing only the sequences that are in OTUs shared by the groups specified. It can only be used with a list and group file not the shared file input.\n";
		helpString += "The output parameter allows you to output the list of names without the group and bin number added. \n";
		helpString += "With this option you can use the names file as an input in get.seqs and remove.seqs commands. To do this enter output=accnos. \n";
		helpString += "The get.sharedseqs command outputs a .names file for each distance level containing a list of sequences in the OTUs shared by the groups specified.\n";
		helpString += "The get.sharedseqs command should be in the following format: get.sharedseqs(list=yourListFile, group=yourGroupFile, label=yourLabels, uniquegroups=yourGroups, fasta=yourFastafile, output=yourOutput).\n";
		helpString += "Example get.sharedseqs(list=amazon.fn.list, label=unique-0.01, group=amazon.groups, uniquegroups=forest-pasture, fasta=amazon.fasta, output=accnos).\n";
		helpString += "The output to the screen is the distance and the number of otus at that distance for the groups you specified.\n";
		helpString += "The default value for label is all labels in your inputfile. The default for groups is all groups in your file.\n";
		helpString += "Note: No spaces between parameter labels (i.e. label), '=' and parameters (i.e.yourLabel).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "GetSharedOTUCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string GetSharedOTUCommand::getOutputPattern(string type) {
    try {
        string pattern = "";

        if (type == "fasta")            {   pattern =  "[filename],[distance],[group],shared.fasta";   }
        else if (type == "accnos")      {   pattern =  "[filename],[distance],[group],accnos";         }
        else if (type == "sharedseqs")  {   pattern =  "[filename],[distance],[group],shared.seqs";    }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "GetSharedOTUCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
GetSharedOTUCommand::GetSharedOTUCommand(){	
	try {
		abort = true; calledHelp = true;
		setParameters();
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
GetSharedOTUCommand::GetSharedOTUCommand(string option)  {
	try {
	
		abort = false; calledHelp = false;   
		unique = true;
		allLines = 1;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
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
			}

			
			//check for required parameters
			listfile = validParameter.validFile(parameters, "list", true);
			if (listfile == "not open") { abort = true; }
			else if (listfile == "not found") { listfile = "";			}
            else {  format = "list"; 	m->setListFile(listfile); }
			
			groupfile = validParameter.validFile(parameters, "group", true);
			if (groupfile == "not open") { abort = true; }	
			else if (groupfile == "not found") { groupfile = ""; }
			else { m->setGroupFile(groupfile); }
            
            sharedfile = validParameter.validFile(parameters, "shared", true);
			if (sharedfile == "not open") { abort = true; }
			else if (sharedfile == "not found") { sharedfile = ""; }
			else { m->setSharedFile(sharedfile); }
            
            fastafile = validParameter.validFile(parameters, "fasta", true);
			if (fastafile == "not open") { abort = true; }
			else if (fastafile == "not found") {  fastafile = "";  }
			else { m->setFastaFile(fastafile); }

            
            if ((sharedfile == "") && (listfile == "")) { //look for currents
                //is there are current file available for either of these?
				//give priority to shared, then list
				sharedfile = m->getSharedFile();
				if (sharedfile != "") {  m->mothurOut("Using " + sharedfile + " as input file for the shared parameter."); m->mothurOutEndLine(); }
				else {
					listfile = m->getListFile();
					if (listfile != "") {  m->mothurOut("Using " + listfile + " as input file for the list parameter."); m->mothurOutEndLine(); }
					else {
						m->mothurOut("No valid current files. You must provide a shared or list file."); m->mothurOutEndLine();
						abort = true;
					}
				}
            }else if ((sharedfile != "") && (listfile != "")) {
                m->mothurOut("You may enter ONLY ONE of the following: shared or list."); m->mothurOutEndLine(); abort = true;
            }
			
            if (listfile != "") {
				if (groupfile == ""){
					groupfile = m->getGroupFile();
					if (groupfile != "") {  m->mothurOut("Using " + groupfile + " as input file for the group parameter."); m->mothurOutEndLine(); }
					else { m->mothurOut("You need to provide a group file if you are going to use the list format."); m->mothurOutEndLine(); abort = true;
					}
				}
			}

			if ((sharedfile != "") && (fastafile != "")) { m->mothurOut("You cannot use the fasta file with the shared file."); m->mothurOutEndLine(); abort = true; }
            
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
			else if (output == "default") { output = ""; }
			
			groups = validParameter.validFile(parameters, "uniquegroups", false);			
			if (groups == "not found") { groups = ""; }
			else { 
				userGroups = "unique." + groups;
				m->splitAtDash(groups, Groups);
                if (Groups.size() > 4) {  userGroups = "unique.selected_groups"; } //if too many groups then the filename becomes too big.
			}
			
			groups = validParameter.validFile(parameters, "sharedgroups", false);			
			if (groups == "not found") { groups = "";  }
			else { 
				userGroups = groups;
				m->splitAtDash(groups, Groups);
                if (Groups.size() > 4) {  userGroups = "selected_groups"; } //if too many groups then the filename becomes too big.
				unique = false;
			}
			
        }

	}
	catch(exception& e) {
		m->errorOut(e, "GetSharedOTUCommand", "GetSharedOTUCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int GetSharedOTUCommand::execute(){
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
        if ( sharedfile != "") { runShared(); }
        else {
            m->setGroups(Groups);
            groupMap = new GroupMap(groupfile);
            int error = groupMap->readMap();
            if (error == 1) { delete groupMap; return 0; }
            
            if (m->control_pressed) { delete groupMap; return 0; }
            
            if (Groups.size() == 0) {
                Groups = groupMap->getNamesOfGroups();
                
                //make string for outputfile name
                userGroups = "unique.";
                for(int i = 0; i < Groups.size(); i++) {  userGroups += Groups[i] + "-";  }
                userGroups = userGroups.substr(0, userGroups.length()-1);
                if (Groups.size() > 4) {  userGroups = "unique.selected_groups"; } //if too many groups then the filename becomes too big.
            }else{
                //sanity check for group names
                SharedUtil util;
                vector<string> namesOfGroups = groupMap->getNamesOfGroups(); 
                util.setGroups(Groups, namesOfGroups);
                groupMap->setNamesOfGroups(namesOfGroups);
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
                    for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); }  outputTypes.clear();
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
            m->clearGroups();  
            
            if (lastlist != NULL) {		delete lastlist;	}
            
            if (m->control_pressed) { outputTypes.clear(); for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); }  delete groupMap; return 0; } 
		}
		//set fasta file as new current fastafile
		string current = "";
		itTypes = outputTypes.find("fasta");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setFastaFile(current); }
		}
		
		if (output == "accnos") {
			itTypes = outputTypes.find("accnos");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setAccnosFile(current); }
			}
		}
		
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
        map<string, string> variables;
        variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(listfile));
        variables["[distance]"] = shared->getLabel();
        variables["[group]"] = userGroups;
		if (output != "accnos") { outputFileNames = getOutputFileName("sharedseqs", variables); }
		else { outputFileNames = getOutputFileName("accnos", variables); }
        
		m->openOutputFile(outputFileNames, outNames);
		
		bool wroteSomething = false;
		int num = 0;
				
		//go through each bin, find out if shared
        vector<string> binLabels = shared->getLabels();
		for (int i = 0; i < shared->getNumBins(); i++) {
			if (m->control_pressed) { outNames.close(); m->mothurRemove(outputFileNames); return 0; }
			
			bool uniqueOTU = true;
			
			map<string, int> atLeastOne;
			for (int f = 0; f < Groups.size(); f++) {
				atLeastOne[Groups[f]] = 0;
			}
			
			vector<string> namesOfSeqsInThisBin;
			
			string names = shared->get(i); 
            vector<string> binNames;
            m->splitAtComma(names, binNames);
			for(int j = 0; j < binNames.size(); j++) {
				string name = binNames[j];
				
				//find group
				string seqGroup = groupMap->getGroup(name);
				if (output != "accnos") {
					namesOfSeqsInThisBin.push_back((name + "|" + seqGroup + "|" + binLabels[i]));
				}else {  namesOfSeqsInThisBin.push_back(name);	}
				
				if (seqGroup == "not found") { m->mothurOut(name + " is not in your groupfile. Please correct."); m->mothurOutEndLine(); exit(1);  }
				
				//is this seq in one of hte groups we care about
				it = groupFinder.find(seqGroup);
				if (it == groupFinder.end()) {  uniqueOTU = false;  } //you have a sequence from a group you don't want
				else {  atLeastOne[seqGroup]++;  }
			}
			
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
			m->mothurRemove(outputFileNames);
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
            variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(fastafile));
			string outputFileFasta = getOutputFileName("fasta", variables);
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
/***********************************************************/
int GetSharedOTUCommand::runShared() {
	try {
        InputData input(sharedfile, "sharedfile");
		vector<SharedRAbundVector*> lookup = input.getSharedRAbundVectors();
		string lastLabel = lookup[0]->getLabel();
        
        if (Groups.size() == 0) {
            Groups = m->getGroups();
            
            //make string for outputfile name
            userGroups = "unique.";
            for(int i = 0; i < Groups.size(); i++) {  userGroups += Groups[i] + "-";  }
            userGroups = userGroups.substr(0, userGroups.length()-1);
            if (Groups.size() > 4) {  userGroups = "unique.selected_groups"; } //if too many groups then the filename becomes too big.
        }else {
            //sanity check for group names
            SharedUtil util;
            vector<string> allGroups = m->getAllGroups();
            util.setGroups(Groups, allGroups);
        }
        
        //put groups in map to find easier
        for(int i = 0; i < Groups.size(); i++) {
            groupFinder[Groups[i]] = Groups[i];
        }

		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;
        
		//as long as you are not at the end of the file or done wih the lines you want
		while((lookup[0] != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			if (m->control_pressed) {
                outputTypes.clear(); for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); }
                for (int i = 0; i < lookup.size(); i++) { delete lookup[i]; } m->clearGroups(); return 0;
			}
            
            
			if(allLines == 1 || labels.count(lookup[0]->getLabel()) == 1){
				m->mothurOut(lookup[0]->getLabel()); 
				process(lookup);
				
				processedLabels.insert(lookup[0]->getLabel());
				userLabels.erase(lookup[0]->getLabel());
			}
			
			if ((m->anyLabelsToProcess(lookup[0]->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
                string saveLabel = lookup[0]->getLabel();
                
                for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }
                lookup = input.getSharedRAbundVectors(lastLabel);
                
                m->mothurOut(lookup[0]->getLabel()); 
                process(lookup);
                
                processedLabels.insert(lookup[0]->getLabel());
                userLabels.erase(lookup[0]->getLabel());
                
                //restore real lastlabel to save below
                lookup[0]->setLabel(saveLabel);
			}
			
			lastLabel = lookup[0]->getLabel();
            
			//get next line to process
			//prevent memory leak
			for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }
			lookup = input.getSharedRAbundVectors();
		}
		
		if (m->control_pressed) {
            outputTypes.clear(); for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); }
            for (int i = 0; i < lookup.size(); i++) { delete lookup[i]; } m->clearGroups(); return 0;
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
            for (int i = 0; i < lookup.size(); i++) {  if (lookup[i] != NULL) {	delete lookup[i];	} }
            lookup = input.getSharedRAbundVectors(lastLabel);
            
            m->mothurOut(lookup[0]->getLabel()); 
            process(lookup);
            for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } 
		}
        
		//reset groups parameter
		m->clearGroups();  
		
		return 0;

    }
	catch(exception& e) {
		m->errorOut(e, "GetSharedOTUCommand", "runShared");
		exit(1);
	}
}
/***********************************************************/
int GetSharedOTUCommand::process(vector<SharedRAbundVector*>& lookup) {
	try {
		
		string outputFileNames;
		if (outputDir == "") { outputDir += m->hasPath(sharedfile); }
        map<string, string> variables;
        variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(sharedfile));
        variables["[distance]"] = lookup[0]->getLabel();
        variables["[group]"] = userGroups;
		if (output != "accnos") { outputFileNames = getOutputFileName("sharedseqs", variables); }
		else { outputFileNames = getOutputFileName("accnos", variables); }
        
        ofstream outNames;
		m->openOutputFile(outputFileNames, outNames);
		
		bool wroteSomething = false;
		int num = 0;
        
		//go through each bin, find out if shared
		for (int i = 0; i < lookup[0]->getNumBins(); i++) {
			if (m->control_pressed) { outNames.close(); m->mothurRemove(outputFileNames); return 0; }
			
			bool uniqueOTU = true;
			map<string, int> atLeastOne;
			for (int f = 0; f < Groups.size(); f++) {  atLeastOne[Groups[f]] = 0;  }
			
			set<string> namesOfGroupsInThisBin;
			
			for(int j = 0; j < lookup.size(); j++) {
				string seqGroup = lookup[j]->getGroup();
                string name = m->currentSharedBinLabels[i];
				
                if (lookup[j]->getAbundance(i) != 0) {
                    if (output != "accnos") {
                        namesOfGroupsInThisBin.insert(name + "|" + seqGroup + "|" + toString(lookup[j]->getAbundance(i)));
                    }else {  namesOfGroupsInThisBin.insert(name);	}
                    
                    //is this seq in one of the groups we care about
                    it = groupFinder.find(seqGroup);
                    if (it == groupFinder.end()) {  uniqueOTU = false;  } //you have sequences from a group you don't want
                    else {  atLeastOne[seqGroup]++;  }
				}
			}
			
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
				for (set<string>::iterator itNames = namesOfGroupsInThisBin.begin(); itNames != namesOfGroupsInThisBin.end(); itNames++) {
					outNames << (*itNames) << endl;
				}
			}
		}
		outNames.close();
		
		if (!wroteSomething) {
			m->mothurRemove(outputFileNames);
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
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "GetSharedOTUCommand", "process");
		exit(1);
	}
}

//**********************************************************************************************************************
