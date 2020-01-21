/*
 *  removerarecommand.cpp
 *  mothur
 *
 *  Created by westcott on 1/21/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "removerarecommand.h"
#include "sequence.hpp"
#include "groupmap.h"

#include "inputdata.h"

//**********************************************************************************************************************
vector<string> RemoveRareCommand::setParameters(){	
	try {
		CommandParameter plist("list", "InputTypes", "", "", "none", "atleast", "none","list",false,false,true); parameters.push_back(plist);
		CommandParameter prabund("rabund", "InputTypes", "", "", "none", "atleast", "none","rabund",false,false,true); parameters.push_back(prabund);
		CommandParameter psabund("sabund", "InputTypes", "", "", "none", "atleast", "none","sabund",false,false,true); parameters.push_back(psabund);
		CommandParameter pshared("shared", "InputTypes", "", "", "none", "atleast", "none","shared",false,false,true); parameters.push_back(pshared);
        CommandParameter pcount("count", "InputTypes", "", "", "CountGroup", "none", "none","count",false,false); parameters.push_back(pcount);
		CommandParameter pgroup("group", "InputTypes", "", "", "CountGroup", "none", "none","group",false,false); parameters.push_back(pgroup);
		CommandParameter pgroups("groups", "String", "", "", "", "", "","",false,false); parameters.push_back(pgroups);
		CommandParameter plabel("label", "String", "", "", "", "", "","",false,false); parameters.push_back(plabel);
		CommandParameter pnseqs("nseqs", "Number", "", "0", "", "", "","",false,true,true); parameters.push_back(pnseqs);
		CommandParameter pbygroup("bygroup", "Boolean", "", "f", "", "", "","",false,false); parameters.push_back(pbygroup);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveRareCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string RemoveRareCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The remove.rare command parameters are list, rabund, sabund, shared, group, count, label, groups, bygroup and nseqs.\n";
		helpString += "The remove.rare command reads one of the following file types: list, rabund, sabund or shared file. It outputs a new file after removing the rare otus.\n";
		helpString += "The groups parameter allows you to specify which of the groups you would like analyzed.  Default=all. You may separate group names with dashes.\n";
		helpString += "The label parameter is used to analyze specific labels in your input. default=all. You may separate label names with dashes.\n";
		helpString += "The bygroup parameter is only valid with the shared file. default=f, meaning remove any OTU that has nseqs or fewer sequences across all groups.\n";
		helpString += "bygroups=T means remove any OTU that has nseqs or fewer sequences in each group (if groupA has 1 sequence and group B has 100 sequences in OTUZ and nseqs=1, then set the groupA count for OTUZ to 0 and keep groupB's count at 100.) \n";
		helpString += "The nseqs parameter allows you to set the cutoff for an otu to be deemed rare. It is required.\n";
		helpString += "The remove.rare command should be in the following format: remove.rare(shared=yourSharedFile, nseqs=yourRareCutoff).\n";
		helpString += "Example remove.rare(shared=amazon.fn.shared, nseqs=2).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveRareCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string RemoveRareCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "rabund")            {   pattern = "[filename],pick,[extension]";    }
        else if (type == "sabund")    {   pattern = "[filename],pick,[extension]";    }
        else if (type == "group")       {   pattern = "[filename],pick,[extension]";    }
        else if (type == "count")       {   pattern = "[filename],pick,[extension]";    }
        else if (type == "list")        {   pattern = "[filename],[tag],pick,[extension]";    }
        else if (type == "shared")      {   pattern = "[filename],[tag],pick,[extension]";    }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "RemoveRareCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
RemoveRareCommand::RemoveRareCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["rabund"] = tempOutNames;
		outputTypes["sabund"] = tempOutNames;
		outputTypes["list"] = tempOutNames;
		outputTypes["group"] = tempOutNames;
        outputTypes["count"] = tempOutNames;
		outputTypes["shared"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveRareCommand", "RemoveRareCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
RemoveRareCommand::RemoveRareCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		allLines = true;
		
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
				if (!validParameter.isValidParameter(it->first, myArray, it->second)) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["rabund"] = tempOutNames;
			outputTypes["sabund"] = tempOutNames;
			outputTypes["list"] = tempOutNames;
			outputTypes["group"] = tempOutNames;
			outputTypes["shared"] = tempOutNames;	
            outputTypes["count"] = tempOutNames;
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){	outputDir = "";		}
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.valid(parameters, "inputdir");		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("list");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["list"] = inputDir + it->second;		}
				}
				
				it = parameters.find("group");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["group"] = inputDir + it->second;		}
				}
				
				it = parameters.find("sabund");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["sabund"] = inputDir + it->second;		}
				}
				
				it = parameters.find("rabund");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["rabund"] = inputDir + it->second;		}
				}
				
				it = parameters.find("shared");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["shared"] = inputDir + it->second;		}
				}
                
                it = parameters.find("count");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["count"] = inputDir + it->second;		}
				}
			}
			
			
			//check for file parameters
			listfile = validParameter.validFile(parameters, "list");
			if (listfile == "not open") { abort = true; }
			else if (listfile == "not found") {  listfile = "";  }
			else { current->setListFile(listfile); }
			
			sabundfile = validParameter.validFile(parameters, "sabund");
			if (sabundfile == "not open") { abort = true; }
			else if (sabundfile == "not found") {  sabundfile = "";  }	
			else { current->setSabundFile(sabundfile); }
			
			rabundfile = validParameter.validFile(parameters, "rabund");
			if (rabundfile == "not open") { abort = true; }
			else if (rabundfile == "not found") {  rabundfile = "";  }				
			else { current->setRabundFile(rabundfile); }
			
			groupfile = validParameter.validFile(parameters, "group");
			if (groupfile == "not open") { groupfile = ""; abort = true; }
			else if (groupfile == "not found") {  groupfile = "";  }	
			else { current->setGroupFile(groupfile); }
			
			sharedfile = validParameter.validFile(parameters, "shared");
			if (sharedfile == "not open") { sharedfile = "";  abort = true; }
			else if (sharedfile == "not found") {  sharedfile = "";  }
			else { current->setSharedFile(sharedfile); }
            
            countfile = validParameter.validFile(parameters, "count");
			if (countfile == "not open") { countfile = ""; abort = true; }
			else if (countfile == "not found") { countfile = "";  }	
			else { current->setCountFile(countfile); }
            			
            if ((groupfile != "") && (countfile != "")) {
                m->mothurOut("[ERROR]: you may only use one of the following: group or count."); m->mothurOutEndLine(); abort=true;
            }
			
			if ((sharedfile == "") && (listfile == "") && (rabundfile == "") && (sabundfile == "")) { 
				//is there are current file available for any of these?
				//give priority to shared, then list, then rabund, then sabund
				//if there is a current shared file, use it
				sharedfile = current->getSharedFile(); 
				if (sharedfile != "") {  m->mothurOut("Using " + sharedfile + " as input file for the shared parameter."); m->mothurOutEndLine(); }
				else { 
					listfile = current->getListFile(); 
					if (listfile != "") {  m->mothurOut("Using " + listfile + " as input file for the list parameter."); m->mothurOutEndLine(); }
					else { 
						rabundfile = current->getRabundFile(); 
						if (rabundfile != "") {  m->mothurOut("Using " + rabundfile + " as input file for the rabund parameter."); m->mothurOutEndLine(); }
						else { 
							sabundfile = current->getSabundFile(); 
							if (sabundfile != "") {  m->mothurOut("Using " + sabundfile + " as input file for the sabund parameter."); m->mothurOutEndLine(); }
							else { 
								m->mothurOut("No valid current files. You must provide a list, sabund, rabund or shared file."); m->mothurOutEndLine(); 
								abort = true;
							}
						}
					}
				}
			} 
			
			groups = validParameter.valid(parameters, "groups");			
			if (groups == "not found") { groups = "all"; }
			util.splitAtDash(groups, Groups);
            if (Groups.size() != 0) { if (Groups[0]== "all") { Groups.clear(); } }
			
			label = validParameter.valid(parameters, "label");			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  util.splitAtDash(label, labels);  allLines = false;  }
				else { allLines = true;  }
			}
			
			string temp = validParameter.valid(parameters, "nseqs");
			if (temp == "not found") { m->mothurOut("nseqs is a required parameter."); m->mothurOutEndLine(); abort = true; }
			else { util.mothurConvert(temp, nseqs); }
			
			temp = validParameter.valid(parameters, "bygroup");	 if (temp == "not found") { temp = "f"; }
			byGroup = util.isTrue(temp);
			
			if (byGroup && (sharedfile == "")) { m->mothurOut("The byGroup parameter is only valid with a shared file.\n"); }
			
			if (((groupfile != "") || (countfile != "")) && (listfile == "")) { m->mothurOut("A group or count file is only valid with a list file.\n");  groupfile = ""; countfile = ""; }
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveRareCommand", "RemoveRareCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int RemoveRareCommand::execute(){
	try {
		
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
		if (m->getControl_pressed()) { return 0; }
		
		//read through the correct file and output lines you want to keep
		if (sabundfile != "")		{		processSabund();	}
		if (rabundfile != "")		{		processRabund();	}
		if (listfile != "")			{		processList();		}
		if (sharedfile != "")		{		processShared();	}
		
		if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); } return 0; }
			
		if (outputNames.size() != 0) {
			m->mothurOutEndLine();
			m->mothurOut("Output File Names: "); m->mothurOutEndLine();
			for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
			m->mothurOutEndLine();
			
			//set rabund file as new current rabundfile
			string currentName = "";
			itTypes = outputTypes.find("rabund");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setRabundFile(currentName); }
			}
			
			itTypes = outputTypes.find("sabund");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setSabundFile(currentName); }
			}
			
			itTypes = outputTypes.find("group");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setGroupFile(currentName); }
			}
			
			itTypes = outputTypes.find("list");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setListFile(currentName); }
			}
			
			itTypes = outputTypes.find("shared");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setSharedFile(currentName); }
			}
            
            itTypes = outputTypes.find("count");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setCountFile(currentName); }
			}
		}
		
		return 0;		
	}
	
	catch(exception& e) {
		m->errorOut(e, "RemoveRareCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************
int RemoveRareCommand::processList(){
	try {
				
		//you must provide a label because the names in the listfile need to be consistent
		string thisLabel = "";
		if (allLines) { m->mothurOut("For the listfile you must select one label, using first label in your listfile."); m->mothurOutEndLine(); }
		else if (labels.size() > 1) { m->mothurOut("For the listfile you must select one label, using " + (*labels.begin()) + "."); m->mothurOutEndLine(); thisLabel = *labels.begin(); }
		else { thisLabel = *labels.begin(); }
		
		InputData input(listfile, "list", nullVector);
		ListVector* list = input.getListVector();
		
		//get first one or the one we want
		if (thisLabel != "") { 	
			//use smart distancing
			set<string> userLabels; userLabels.insert(thisLabel);
			set<string> processedLabels;
			string lastLabel = list->getLabel();
			while((list != NULL) && (userLabels.size() != 0)) {
				if(userLabels.count(list->getLabel()) == 1){
					processedLabels.insert(list->getLabel());
					userLabels.erase(list->getLabel());
					break;
				}
				
				if ((util.anyLabelsToProcess(list->getLabel(), userLabels, "") ) && (processedLabels.count(lastLabel) != 1)) {
					processedLabels.insert(list->getLabel());
					userLabels.erase(list->getLabel());
					delete list;
					list = input.getListVector(lastLabel);
					break;
				}
				lastLabel = list->getLabel();
				delete list;
				list = input.getListVector();
			}
			if (userLabels.size() != 0) { 
				m->mothurOut("Your file does not include the label " + thisLabel + ". I will use " + lastLabel + ".");  m->mothurOutEndLine();
				list = input.getListVector(lastLabel); 
			}
		}
        
        string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += util.hasPath(listfile);  }
        map<string, string> variables;
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(listfile));
        variables["[extension]"] = util.getExtension(listfile);
        variables["[tag]"] = list->getLabel();
		string outputFileName = getOutputFileName("list", variables);
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(groupfile));
        variables["[extension]"] = util.getExtension(groupfile);
		string outputGroupFileName = getOutputFileName("group", variables);
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(countfile));
        variables["[extension]"] = util.getExtension(countfile);
        string outputCountFileName = getOutputFileName("count", variables);
        
		ofstream out, outGroup;
		util.openOutputFile(outputFileName, out);
		
		bool wroteSomething = false;

		
		//if groupfile is given then use it
		GroupMap* groupMap;
        CountTable ct;
        CountTable newCountTable; //instead of removing rare, fill new count table with "good" seqs
        bool selectedGroups = true;
		if (groupfile != "") { 
			groupMap = new GroupMap(groupfile); groupMap->readMap(); 
            if (Groups.size() == 0) { Groups = groupMap->getNamesOfGroups(); selectedGroups = false; }
			util.openOutputFile(outputGroupFileName, outGroup);
		}else if (countfile != "") {
            ct.readTable(countfile, true, false);
            if (ct.hasGroupInfo()) {
                vector<string> namesGroups = ct.getNamesOfGroups();
                if (Groups.size() == 0) { Groups = ct.getNamesOfGroups(); selectedGroups = false; }
                for (int i = 0; i < namesGroups.size(); i++) { newCountTable.addGroup(namesGroups[i]); }
            }
        }
        
		if (list != NULL) {
            
            vector<string> binLabels = list->getLabels();
            vector<string> newLabels;
            
			//make a new list vector
			ListVector newList;
			newList.setLabel(list->getLabel());

			//for each bin
			for (int i = 0; i < list->getNumBins(); i++) {
				if (m->getControl_pressed()) {  if (groupfile != "") { delete groupMap; outGroup.close(); util.mothurRemove(outputGroupFileName); } out.close();  util.mothurRemove(outputFileName);  return 0; }
                
				//parse out names that are in accnos file
				string binnames = list->get(i);
				vector<string> names;
                vector<string> newNames;
				util.splitAtComma(binnames, names);
                int binsize = names.size();
				
				vector<string> newGroupFile;
				if (groupfile != "") {
					for(int k = 0; k < names.size(); k++) {
						string group = groupMap->getGroup(names[k]);
						
                        if (selectedGroups) {
                            if (util.inUsersGroups(group, Groups)) {
                                newGroupFile.push_back(names[k] + "\t" + group);
                                newNames.push_back(names[k]);
                            }
                        }else {
                            newGroupFile.push_back(names[k] + "\t" + group);
                            newNames.push_back(names[k]);
                        }
					}
					names = newNames; binsize = names.size();
				}else if (countfile != "") {
                    binsize = 0;
					for(int k = 0; k < names.size(); k++) {
                        if (ct.hasGroupInfo()) {
                            if (selectedGroups) {
                                vector<string> thisSeqsGroups = ct.getGroups(names[k]);
                                vector<int> thisGroupCounts = ct.getGroupCounts(names[k]);
                            
                                int thisSeqsCount = 0;
                                for (int n = 0; n < thisSeqsGroups.size(); n++) {
                                    if (util.inUsersGroups(thisSeqsGroups[n], Groups)) {
                                        thisSeqsCount += thisGroupCounts[n];
                                    }
                                }
                                binsize += thisSeqsCount;
                                //if you don't have any seqs from the groups the user wants, then remove you.
                                if (thisSeqsCount == 0) { newGroupFile.push_back(names[k]); }
                                else { newNames.push_back(names[k]); }
                            }else { //all groups
                                binsize += ct.getNumSeqs(names[k]);
                                newNames.push_back(names[k]);
                            }
                        }else {
                            binsize += ct.getNumSeqs(names[k]); 
                            newNames.push_back(names[k]);
                        }
					}
                }

				if (binsize > nseqs) { //keep bin
                    string saveBinNames = util.getStringFromVector(newNames, ",");
					newList.push_back(saveBinNames);
                    newLabels.push_back(binLabels[i]);
					if (groupfile != "") {  for(int k = 0; k < newGroupFile.size(); k++) { outGroup << newGroupFile[k] << endl; }  }
                    else if (countfile != "") {
                        for(int k = 0; k < newNames.size(); k++) {
                            if (ct.hasGroupInfo()) {
                                vector<int> groupCounts = ct.getGroupCounts(newNames[k]);
                                newCountTable.push_back(newNames[k], groupCounts);
                            }else {
                                int reps = ct.getNumSeqs(newNames[k]);
                                newCountTable.push_back(newNames[k], reps);
                            }
                        }
                    }
                }
			}
			
			//print new listvector
			if (newList.getNumBins() != 0) {
				wroteSomething = true;
				newList.setLabels(newLabels);
                newList.print(out, false);
			}
		}	
		
		out.close();
		if (groupfile != "") { outGroup.close(); outputTypes["group"].push_back(outputGroupFileName); outputNames.push_back(outputGroupFileName); }
        if (countfile != "") { 
            if (newCountTable.hasGroupInfo()) {
                vector<string> allGroups = newCountTable.getNamesOfGroups();
                for (int i = 0; i < allGroups.size(); i++) {
                    if (!util.inUsersGroups(allGroups[i], Groups)) { newCountTable.removeGroup(allGroups[i]); }
                }
            }
            newCountTable.printTable(outputCountFileName, true);
            outputTypes["count"].push_back(outputCountFileName); outputNames.push_back(outputCountFileName); 
        }
		
		if (wroteSomething == false) {  m->mothurOut("Your file contains only rare sequences."); m->mothurOutEndLine();  }
		outputTypes["list"].push_back(outputFileName); outputNames.push_back(outputFileName);
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveRareCommand", "processList");
		exit(1);
	}
}
//**********************************************************************************************************************
int RemoveRareCommand::processSabund(){
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += util.hasPath(sabundfile);  }
        map<string, string> variables; 
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(sabundfile));
        variables["[extension]"] = util.getExtension(sabundfile);
		string outputFileName = getOutputFileName("sabund", variables);
		outputTypes["sabund"].push_back(outputFileName); outputNames.push_back(outputFileName);

		ofstream out;
		util.openOutputFile(outputFileName, out);
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		InputData input(sabundfile, "sabund", nullVector);
		SAbundVector* sabund = input.getSAbundVector();
		string lastLabel = sabund->getLabel();
		set<string> processedLabels;
		set<string> userLabels = labels;
		
		while((sabund != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			
			if (m->getControl_pressed()) { delete sabund; out.close(); return 0; }
			
			if(allLines == 1 || labels.count(sabund->getLabel()) == 1){			
				
				m->mothurOut(sabund->getLabel()); m->mothurOutEndLine();
				processedLabels.insert(sabund->getLabel());
				userLabels.erase(sabund->getLabel());
				
				if (sabund->getMaxRank() > nseqs) {
					for(int i = 1; i <=nseqs; i++) {  sabund->set(i, 0); }
				}else {	sabund->clear(); }
				
				if (sabund->getNumBins() > 0) { sabund->print(out); }
			}
			
			if ((util.anyLabelsToProcess(sabund->getLabel(), userLabels, "") ) && (processedLabels.count(lastLabel) != 1)) {
				string saveLabel = sabund->getLabel();
				
				delete sabund;
				sabund = input.getSAbundVector(lastLabel);
				
				m->mothurOut(sabund->getLabel()); m->mothurOutEndLine();
				processedLabels.insert(sabund->getLabel());
				userLabels.erase(sabund->getLabel());
				
				if (sabund->getMaxRank() > nseqs) {
					for(int i = 1; i <=nseqs; i++) {  sabund->set(i, 0); }
				}else {	sabund->clear(); }
				
				if (sabund->getNumBins() > 0) { sabund->print(out); }
								
				//restore real lastlabel to save below
				sabund->setLabel(saveLabel);
			}		
			
			lastLabel = sabund->getLabel();			
			
			delete sabund;
			sabund = input.getSAbundVector();
		}
		
		if (m->getControl_pressed()) {  out.close(); return 0; }	
		
		//output error messages about any remaining user labels
		bool needToRun = false;
		for (set<string>::iterator it = userLabels.begin(); it != userLabels.end(); it++) {
			m->mothurOut("Your file does not include the label " + *it); 
            if (processedLabels.count(lastLabel) != 1)  { m->mothurOut(". I will use " + lastLabel + ".\n"); needToRun = true;  }
			else                                        { m->mothurOut(". Please refer to " + lastLabel + ".\n");               }
		}
		
		//run last label if you need to
		if (needToRun )  {
			if (sabund != NULL) {	delete sabund;	}
			sabund = input.getSAbundVector(lastLabel);
			
			m->mothurOut(sabund->getLabel()); m->mothurOutEndLine();
			
			if (sabund->getMaxRank() > nseqs) {
				for(int i = 1; i <=nseqs; i++) {  sabund->set(i, 0); }
			}else {	sabund->clear(); }
			
			if (sabund->getNumBins() > 0) { sabund->print(out); }
			
			delete sabund;
		}
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveRareCommand", "processSabund");
		exit(1);
	}
}
//**********************************************************************************************************************
int RemoveRareCommand::processRabund(){
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += util.hasPath(rabundfile);  }
        map<string, string> variables; 
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(rabundfile));
        variables["[extension]"] = util.getExtension(rabundfile);
		string outputFileName = getOutputFileName("rabund", variables);
		outputTypes["rabund"].push_back(outputFileName); outputNames.push_back(outputFileName);
		
		ofstream out;
		util.openOutputFile(outputFileName, out);
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		InputData input(rabundfile, "rabund", nullVector);
		RAbundVector* rabund = input.getRAbundVector();
		string lastLabel = rabund->getLabel();
		set<string> processedLabels;
		set<string> userLabels = labels;
		
		while((rabund != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			
			if (m->getControl_pressed()) { delete rabund; out.close(); return 0; }
			
			if(allLines == 1 || labels.count(rabund->getLabel()) == 1){			
				
				m->mothurOut(rabund->getLabel()); m->mothurOutEndLine();
				processedLabels.insert(rabund->getLabel());
				userLabels.erase(rabund->getLabel());
				
				RAbundVector newRabund; newRabund.setLabel(rabund->getLabel());
				for (int i = 0; i < rabund->getNumBins(); i++) {
					if (rabund->get(i) > nseqs) {
						newRabund.push_back(rabund->get(i));
					}
				}
				if (newRabund.getNumBins() > 0) { newRabund.print(out); }
			}
			
			if ((util.anyLabelsToProcess(rabund->getLabel(), userLabels, "") ) && (processedLabels.count(lastLabel) != 1)) {
				string saveLabel = rabund->getLabel();
				
				delete rabund;
				rabund = input.getRAbundVector(lastLabel);
				
				m->mothurOut(rabund->getLabel()); m->mothurOutEndLine();
				processedLabels.insert(rabund->getLabel());
				userLabels.erase(rabund->getLabel());
				
				RAbundVector newRabund; newRabund.setLabel(rabund->getLabel());
				for (int i = 0; i < rabund->getNumBins(); i++) {
					if (rabund->get(i) > nseqs) {
						newRabund.push_back(rabund->get(i));
					}
				}
				if (newRabund.getNumBins() > 0) { newRabund.print(out); }				
				
				//restore real lastlabel to save below
				rabund->setLabel(saveLabel);
			}		
			
			lastLabel = rabund->getLabel();			
			
			delete rabund;
			rabund = input.getRAbundVector();
		}
		
		if (m->getControl_pressed()) {  out.close(); return 0; }	
		
		//output error messages about any remaining user labels
		bool needToRun = false;
		for (set<string>::iterator it = userLabels.begin(); it != userLabels.end(); it++) {
			m->mothurOut("Your file does not include the label " + *it); 
            if (processedLabels.count(lastLabel) != 1)  { m->mothurOut(". I will use " + lastLabel + ".\n"); needToRun = true;  }
			else                                        { m->mothurOut(". Please refer to " + lastLabel + ".\n");               }
		}
		
		//run last label if you need to
		if (needToRun )  {
			if (rabund != NULL) {	delete rabund;	}
			rabund = input.getRAbundVector(lastLabel);
			
			m->mothurOut(rabund->getLabel()); m->mothurOutEndLine();
			
			RAbundVector newRabund; newRabund.setLabel(rabund->getLabel());
			for (int i = 0; i < rabund->getNumBins(); i++) {
				if (rabund->get(i) > nseqs) {
					newRabund.push_back(rabund->get(i));
				}
			}
			if (newRabund.getNumBins() > 0) { newRabund.print(out); }	
			
			delete rabund;
		}
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveRareCommand", "processRabund");
		exit(1);
	}
}
//**********************************************************************************************************************
int RemoveRareCommand::processShared(){
	try {
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		InputData input(sharedfile, "sharedfile", Groups);
		set<string> processedLabels;
        set<string> userLabels = labels;
        string lastLabel = "";
        
        SharedRAbundVectors* lookup = util.getNextShared(input, allLines, userLabels, processedLabels, lastLabel);
        Groups = lookup->getNamesGroups();
        
        while (lookup != NULL) {
            
            if (m->getControl_pressed()) { delete lookup; break; }
            
            processLookup(lookup); delete lookup;
            
            lookup = util.getNextShared(input, allLines, userLabels, processedLabels, lastLabel);
        }
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveRareCommand", "processShared");
		exit(1);
	}
}
//**********************************************************************************************************************
int RemoveRareCommand::processLookup(SharedRAbundVectors*& lookup){
	try {
		
        string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += util.hasPath(sharedfile);  }
        map<string, string> variables;
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(sharedfile));
        variables["[extension]"] = util.getExtension(sharedfile);
        variables["[tag]"] = lookup->getLabel();
		string outputFileName = getOutputFileName("shared", variables);
		outputTypes["shared"].push_back(outputFileName); outputNames.push_back(outputFileName);
		
		ofstream out;
		util.openOutputFile(outputFileName, out);
        
		vector<SharedRAbundVector> newRabunds;  newRabunds.resize(lookup->size());
        vector<string> headers;
        vector<string> namesOfGroups = lookup->getNamesGroups();
		for (int i = 0; i < namesOfGroups.size(); i++) {
			newRabunds[i].setGroup(namesOfGroups[i]);
			newRabunds[i].setLabel(lookup->getLabel());
		}
		
        vector<SharedRAbundVector*> data = lookup->getSharedRAbundVectors();
        vector<string> currentLabels = lookup->getOTUNames();
		if (byGroup) {
			
			//for each otu
			for (int i = 0; i < lookup->getNumBins(); i++) {
				bool allZero = true;
				
				if (m->getControl_pressed()) { out.close(); return 0; }
				
				//for each group
                vector<int> abunds = lookup->getOTU(i);
				for (int j = 0; j < abunds.size(); j++) {
					
					//are you rare?
                    if (abunds[j] > nseqs) { allZero = false; }
                    else { abunds[j] = 0; }
				}
				
				//eliminates zero otus
				if (allZero) { }
                else {
                    for (int j = 0; j < abunds.size(); j++) {  newRabunds[j].push_back(abunds[j]); }
                    headers.push_back(currentLabels[i]);
                }
			}
		}else {
			//for each otu
			for (int i = 0; i < lookup->getNumBins(); i++) {
				
				if (m->getControl_pressed()) { out.close(); return 0; }
				
                int totalAbund = lookup->getOTUTotal(i);
				
				//eliminates otus below rare cutoff
				if (totalAbund <= nseqs) {  } //ignore
                else {
                    headers.push_back(currentLabels[i]);
                    for (int j = 0; j < data.size(); j++) { newRabunds[j].push_back(data[j]->get(i)); }
                }
			}
		}
		
		//do we have any otus above the rare cutoff
		if (newRabunds[0].getNumBins() != 0) {
            out << "label\tGroup\tnumOtus";
            for (int j = 0; j < headers.size(); j++) { out << '\t' << headers[j]; }
            out << endl;
			for (int j = 0; j < newRabunds.size(); j++) { newRabunds[j].print(out);  }
		}
		
        out.close();
        
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveRareCommand", "processLookup");
		exit(1);
	}
}
//**********************************************************************************************************************




