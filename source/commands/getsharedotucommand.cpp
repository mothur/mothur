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
vector<string> GetSharedOTUCommand::setParameters(){	
	try {
		CommandParameter pfasta("fasta", "InputTypes", "", "", "sharedFasta", "none", "none","fasta",false,false); parameters.push_back(pfasta);
		CommandParameter pgroup("group", "InputTypes", "", "", "none", "GroupCount", "groupList","",false,false,true); parameters.push_back(pgroup);
        CommandParameter pcount("count", "InputTypes", "", "", "none", "GroupCount", "none","",false,false); parameters.push_back(pcount);
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
        helpString += "The count parameter allows you to provide a count file containing the group info for the list file.\n";
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
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
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
	
		abort = false; calledHelp = false;   userGroups = "";
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
			outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){	outputDir = "";		}
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.valid(parameters, "inputdir");		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("fasta");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["fasta"] = inputDir + it->second;		}
				}
				
				it = parameters.find("list");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["list"] = inputDir + it->second;		}
				}
                
                it = parameters.find("shared");
				//user has given a template file
				if(it != parameters.end()){
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["shared"] = inputDir + it->second;		}
				}
				
				it = parameters.find("group");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["group"] = inputDir + it->second;		}
				}
                
                it = parameters.find("count");
                //user has given a template file
                if(it != parameters.end()){
                    path = util.hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["count"] = inputDir + it->second;		}
                }
			}

			
			//check for required parameters
			listfile = validParameter.validFile(parameters, "list");
			if (listfile == "not open") { abort = true; }
			else if (listfile == "not found") { listfile = "";			}
            else {  format = "list"; 	current->setListFile(listfile); }
			
			groupfile = validParameter.validFile(parameters, "group");
			if (groupfile == "not open") { abort = true; }	
			else if (groupfile == "not found") { groupfile = ""; }
			else { current->setGroupFile(groupfile); }
            
            sharedfile = validParameter.validFile(parameters, "shared");
			if (sharedfile == "not open") { abort = true; }
			else if (sharedfile == "not found") { sharedfile = ""; }
			else { current->setSharedFile(sharedfile); }
            
            fastafile = validParameter.validFile(parameters, "fasta");
			if (fastafile == "not open") { abort = true; }
			else if (fastafile == "not found") {  fastafile = "";  }
			else { current->setFastaFile(fastafile); }

            countfile = validParameter.validFile(parameters, "count");
            if (countfile == "not open") { countfile = ""; abort = true; }
            else if (countfile == "not found") { countfile = ""; }
            else {
                current->setCountFile(countfile);
                CountTable temp;
                if (!temp.testGroups(countfile)) { m->mothurOut("[ERROR]: Your count file does not have group info, aborting."); m->mothurOutEndLine(); abort=true; }
            }
            
            if ((sharedfile == "") && (listfile == "")) { //look for currents
                //is there are current file available for either of these?
				//give priority to shared, then list
				sharedfile = current->getSharedFile();
				if (sharedfile != "") {  m->mothurOut("Using " + sharedfile + " as input file for the shared parameter."); m->mothurOutEndLine(); }
				else {
					listfile = current->getListFile();
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
                if ((groupfile == "") && (countfile == "")) {
                    groupfile = current->getGroupFile();
                    if (groupfile != "") {  m->mothurOut("Using " + groupfile + " as input file for the group parameter."); m->mothurOutEndLine(); }
                    else {
                        countfile = current->getCountFile();
                        if (countfile != "") {  m->mothurOut("Using " + countfile + " as input file for the count parameter."); m->mothurOutEndLine(); }
                        else {
                            m->mothurOut("You need to provide a groupfile or countfile if you are going to use the list format."); m->mothurOutEndLine();
                            abort = true;
                        }
                    }
                }
			}

			if ((sharedfile != "") && (fastafile != "")) { m->mothurOut("You cannot use the fasta file with the shared file."); m->mothurOutEndLine(); abort = true; }
            
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			label = validParameter.valid(parameters, "label");			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  util.splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
			
			output = validParameter.valid(parameters, "output");
			if (output == "not found") { output = ""; }
			else if (output == "default") { output = ""; }
			
			groups = validParameter.valid(parameters, "uniquegroups");
			if (groups == "not found") { groups = ""; }
			else { 
				userGroups = "unique." + groups;
				util.splitAtDash(groups, Groups);
                if (Groups.size() != 0) { if (Groups[0]== "all") { Groups.clear(); } }
                if (Groups.size() > 4) {  userGroups = "unique.selected_groups"; } //if too many groups then the filename becomes too big.
			}
			
			groups = validParameter.valid(parameters, "sharedgroups");
			if (groups == "not found") { groups = "";  }
			else { 
				userGroups = groups;
				util.splitAtDash(groups, Groups);
                if (Groups.size() != 0) { if (Groups[0]== "all") { Groups.clear(); } }
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
		
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
        if ( sharedfile != "") { runShared(); }
        else {
            if (groupfile != "") {
                groupMap = new GroupMap(groupfile);
                
                int groupError = groupMap->readMap();
                if (groupError == 1) { delete groupMap; return 0; }
                vector<string> allGroups = groupMap->getNamesOfGroups();
            }else{
                ct = new CountTable();
                ct->readTable(countfile, true, false);
            }
            
            if (m->getControl_pressed()) { delete groupMap; return 0; }
            
            if (Groups.size() == 0) {
                if (groupfile != "") { Groups = groupMap->getNamesOfGroups(); }
                else {  Groups = ct->getNamesOfGroups();  }
                
                //make string for outputfile name
                userGroups = "unique.";
                for(int i = 0; i < Groups.size(); i++) {  userGroups += Groups[i] + "-";  }
                userGroups = userGroups.substr(0, userGroups.length()-1);
                if (Groups.size() > 4) {  userGroups = "unique.selected_groups"; } //if too many groups then the filename becomes too big.
            }
        
            //put groups in map to find easier
            for(int i = 0; i < Groups.size(); i++) { groupFinder[Groups[i]] = Groups[i]; }
            
            if (fastafile != "") {
                ifstream inFasta;
                util.openInputFile(fastafile, inFasta);
                
                while(!inFasta.eof()) {
                    if (m->getControl_pressed()) { outputTypes.clear(); inFasta.close(); delete groupMap; return 0; }
                    
                    Sequence seq(inFasta); util.gobble(inFasta);
                    if (seq.getName() != "") {  seqs.push_back(seq);   }
                }
                inFasta.close();
            }
            
            InputData input(listfile, "list", nullVector);
            ListVector* list = input.getListVector();
            string lastLabel = "";
            
            //if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
            set<string> processedLabels;
            set<string> userLabels = labels;
            
            
            while((list != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
                
                if (m->getControl_pressed()) {
                    for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); }  outputTypes.clear();
                    if (groupfile != "") { delete groupMap; }else { delete ct; } return 0;
                }
                
                if(allLines == 1 || labels.count(list->getLabel()) == 1){
                    m->mothurOut(list->getLabel());
                    process(list);
                
                    processedLabels.insert(list->getLabel());
                    userLabels.erase(list->getLabel());
                }
                
                if ((util.anyLabelsToProcess(list->getLabel(), userLabels, "") ) && (processedLabels.count(lastLabel) != 1)) {
                    string saveLabel = list->getLabel();
                    
                    delete list;
                    list = input.getListVector(lastLabel);
                    
                    m->mothurOut(list->getLabel());
                    process(list);
                    
                    processedLabels.insert(list->getLabel());
                    userLabels.erase(list->getLabel());
                    
                    //restore real lastlabel to save below
                    list->setLabel(saveLabel);
                }
                
                lastLabel = list->getLabel();			
                
                delete list;
                list = input.getListVector();
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
            if (needToRun )  {
                if (list != NULL) {		delete list;	}
                list = input.getListVector(lastLabel);
                
                process(list);
                delete list;
            }
            
            if (m->getControl_pressed()) { outputTypes.clear(); for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); }  if (groupfile != "") { delete groupMap; }else { delete ct; } return 0; } 
		}
		//set fasta file as new current fastafile
		string currentName = "";
		itTypes = outputTypes.find("fasta");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setFastaFile(currentName); }
		}
		
		if (output == "accnos") {
			itTypes = outputTypes.find("accnos");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setAccnosFile(currentName); }
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
		
		if (outputDir == "") { outputDir += util.hasPath(listfile); }
        map<string, string> variables;
        variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(listfile));
        variables["[distance]"] = shared->getLabel();
        variables["[group]"] = userGroups;
		if (output != "accnos") { outputFileNames = getOutputFileName("sharedseqs", variables); }
		else { outputFileNames = getOutputFileName("accnos", variables); }
        
		util.openOutputFile(outputFileNames, outNames);
		
		bool wroteSomething = false;
		int num = 0;
				
		//go through each bin, find out if shared
        vector<string> binLabels = shared->getLabels();
		for (int i = 0; i < shared->getNumBins(); i++) {
			if (m->getControl_pressed()) { outNames.close(); util.mothurRemove(outputFileNames); return 0; }
			
			bool uniqueOTU = true;
			
			map<string, int> atLeastOne;
			for (int f = 0; f < Groups.size(); f++) { atLeastOne[Groups[f]] = 0; }
			
			vector<string> namesOfSeqsInThisBin;
			
			string names = shared->get(i); 
            vector<string> binNames;
            util.splitAtComma(names, binNames);
			for(int j = 0; j < binNames.size(); j++) {
				string name = binNames[j];
				
				//find group
                string seqGroup = "not found"; vector<string> seqsGroups;
                if (groupfile != "") {  seqGroup = groupMap->getGroup(name); }
                else {
                    seqsGroups = ct->getGroups(name);
                    seqGroup = util.getStringFromVector(seqsGroups, "-");
                }
                
				if (output != "accnos") {
					namesOfSeqsInThisBin.push_back((name + "|" + seqGroup + "|" + binLabels[i]));
				}else {  namesOfSeqsInThisBin.push_back(name);	}
				
				if (seqGroup == "not found") { m->mothurOut(name + " is not in your groupfile. Please correct."); m->mothurOutEndLine(); exit(1);  }
				
                if (groupfile != "") {
                    //is this seq in one of hte groups we care about
                    it = groupFinder.find(seqGroup);
                    if (it == groupFinder.end()) {  uniqueOTU = false;  } //you have a sequence from a group you don't want
                    else {  atLeastOne[seqGroup]++;  }
                }else {
                    for (int k = 0; k < seqsGroups.size(); k++) {
                        //is this seq in one of hte groups we care about
                        it = groupFinder.find(seqsGroups[k]);
                        if (it == groupFinder.end()) {  uniqueOTU = false;  } //you have a sequence from a group you don't want
                        else {  atLeastOne[seqsGroups[k]]++;  }
                    }
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
			util.mothurRemove(outputFileNames);
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
			if (outputDir == "") { outputDir += util.hasPath(fastafile); }
            variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(fastafile));
			string outputFileFasta = getOutputFileName("fasta", variables);
			ofstream outFasta;
			util.openOutputFile(outputFileFasta, outFasta);
			outputNames.push_back(outputFileFasta); outputTypes["fasta"].push_back(outputFileFasta);
			
			for (int k = 0; k < seqs.size(); k++) {
				if (m->getControl_pressed()) { outFasta.close(); return 0; }
			
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
        InputData input(sharedfile, "sharedfile", Groups);
		SharedRAbundVectors* lookup = input.getSharedRAbundVectors();
        Groups = lookup->getNamesGroups();
		string lastLabel = lookup->getLabel();
        
        if (userGroups == "") {
            //make string for outputfile name
            userGroups = "unique.";
            for(int i = 0; i < Groups.size(); i++) {  userGroups += Groups[i] + "-";  }
            userGroups = userGroups.substr(0, userGroups.length()-1);
            if (Groups.size() > 4) {  userGroups = "unique.selected_groups"; } //if too many groups then the filename becomes too big.
        }
        
        //put groups in map to find easier
        for(int i = 0; i < Groups.size(); i++) { groupFinder[Groups[i]] = Groups[i];}

		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;
        
		//as long as you are not at the end of the file or done wih the lines you want
		while((lookup != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			if (m->getControl_pressed()) {
                outputTypes.clear(); for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); }
                delete lookup;  return 0;
			}
            
            
			if(allLines == 1 || labels.count(lookup->getLabel()) == 1){
				m->mothurOut(lookup->getLabel());
				process(lookup);
				
				processedLabels.insert(lookup->getLabel());
				userLabels.erase(lookup->getLabel());
			}
			
			if ((util.anyLabelsToProcess(lookup->getLabel(), userLabels, "") ) && (processedLabels.count(lastLabel) != 1)) {
                string saveLabel = lookup->getLabel();
                
                delete lookup;
                lookup = input.getSharedRAbundVectors(lastLabel);
                
                m->mothurOut(lookup->getLabel());
                process(lookup);
                
                processedLabels.insert(lookup->getLabel());
                userLabels.erase(lookup->getLabel());
                
                //restore real lastlabel to save below
                lookup->setLabels(saveLabel);
			}
			
			lastLabel = lookup->getLabel();
            
			//get next line to process
			//prevent memory leak
			delete lookup;
			lookup = input.getSharedRAbundVectors();
		}
		
		if (m->getControl_pressed()) {
            outputTypes.clear(); for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); }
            delete lookup;  return 0;
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
		if (needToRun )  {
            delete lookup;
            lookup = input.getSharedRAbundVectors(lastLabel);
            
            m->mothurOut(lookup->getLabel());
            process(lookup);
            delete lookup;
		}
        
		//reset groups parameter
		  
		
		return 0;

    }
	catch(exception& e) {
		m->errorOut(e, "GetSharedOTUCommand", "runShared");
		exit(1);
	}
}
/***********************************************************/
int GetSharedOTUCommand::process(SharedRAbundVectors*& lookup) {
	try {
		
		string outputFileNames;
		if (outputDir == "") { outputDir += util.hasPath(sharedfile); }
        map<string, string> variables;
        variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(sharedfile));
        variables["[distance]"] = lookup->getLabel();
        variables["[group]"] = userGroups;
		if (output != "accnos") { outputFileNames = getOutputFileName("sharedseqs", variables); }
		else { outputFileNames = getOutputFileName("accnos", variables); }
        
        ofstream outNames;
		util.openOutputFile(outputFileNames, outNames);
		
		bool wroteSomething = false;
		int num = 0;
        
		//go through each bin, find out if shared
		for (int i = 0; i < lookup->getNumBins(); i++) {
			if (m->getControl_pressed()) { outNames.close(); util.mothurRemove(outputFileNames); return 0; }
			
			bool uniqueOTU = true;
			map<string, int> atLeastOne;
			for (int f = 0; f < Groups.size(); f++) {  atLeastOne[Groups[f]] = 0;  }
			
			set<string> namesOfGroupsInThisBin;
			
            vector<string> groupNames = lookup->getNamesGroups();
			for(int j = 0; j < lookup->size(); j++) {
				string seqGroup = groupNames[j];
                string name = lookup->getOTUName(i);
                int abund = lookup->get(i, seqGroup);
				
                if (abund != 0) {
                    if (output != "accnos") {
                        namesOfGroupsInThisBin.insert(name + "|" + seqGroup + "|" + toString(abund));
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
			util.mothurRemove(outputFileNames);
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
