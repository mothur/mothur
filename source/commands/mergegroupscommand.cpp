/*
 *  mergegroupscommand.cpp
 *  mothur
 *
 *  Created by westcott on 1/24/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "mergegroupscommand.h"
#include "sharedutilities.h"
#include "counttable.h"
#include "removeseqscommand.h"

//**********************************************************************************************************************
vector<string> MergeGroupsCommand::setParameters(){	
	try {
		CommandParameter pshared("shared", "InputTypes", "", "", "none", "sharedGroup", "none","shared",false,false,true); parameters.push_back(pshared);
		CommandParameter pgroup("group", "InputTypes", "", "", "CountGroup", "sharedGroup", "none","group",false,false,true); parameters.push_back(pgroup);
        CommandParameter pcount("count", "InputTypes", "", "", "CountGroup", "sharedGroup", "countfasta","count",false,false,true); parameters.push_back(pcount);
		CommandParameter pdesign("design", "InputTypes", "", "", "none", "none", "none","",false,true,true); parameters.push_back(pdesign);
        CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "countfasta","fasta",false,false,true); parameters.push_back(pfasta);
        CommandParameter pmethod("method", "Multiple", "sum-average-median", "sum", "", "", "","",false,false, true); parameters.push_back(pmethod);
		CommandParameter plabel("label", "String", "", "", "", "", "","",false,false); parameters.push_back(plabel);
		CommandParameter pgroups("groups", "String", "", "", "", "", "","",false,false); parameters.push_back(pgroups);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "MergeGroupsCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string MergeGroupsCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The merge.groups command input files are shared, group, count, fasta and a design file.  It reads the design file and merges the groups in the other files accordingly.\n";
		helpString += "The design parameter allows you to assign your groups to sets. It is required. \n";
        helpString += "The fasta parameter allows you to provide a fasta file associated with your count file.  This is used if you are using the median method, so that sequences that are entirely removed from the counttable will also be removed from the fasta file. \n";
		helpString += "The groups parameter allows you to specify which of the groups in your shared or group file you would like included. The group names are separated by dashes. By default all groups are selected.\n";
		helpString += "The label parameter allows you to select what distance levels you would like, and are also separated by dashes.\n";
        helpString += "The groups parameter allows you to select groups you would like, and are also separated by dashes.\n";
        helpString += "The method parameter allows you to select method you would like to use to merge the groups. Options are sum, average and median. Default=sum.\n";
		helpString += "The merge.groups command should be in the following format: merge.groups(design=yourDesignFile, shared=yourSharedFile).\n";
		helpString += "Example merge.groups(design=temp.design, groups=A-B-C, shared=temp.shared).\n";
		helpString += "The default value for groups is all the groups in your sharedfile, and all labels in your inputfile will be used.\n";
		helpString += "Note: No spaces between parameter labels (i.e. groups), '=' and parameters (i.e.yourGroups).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "MergeGroupsCommand", "getHelpString");
		exit(1);
	}
}

//**********************************************************************************************************************
string MergeGroupsCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "shared") {  pattern = "[filename],merge,[extension]"; } 
        else if (type == "group") {  pattern = "[filename],merge,[extension]"; }
        else if (type == "count") {  pattern = "[filename],merge,[extension]"; }
        else if (type == "fasta") {  pattern = "[filename],merge,[extension]"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "MergeGroupsCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
MergeGroupsCommand::MergeGroupsCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["shared"] = tempOutNames;
		outputTypes["group"] = tempOutNames;
        outputTypes["count"] = tempOutNames;
        outputTypes["fasta"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "MergeGroupsCommand", "MergeGroupsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

MergeGroupsCommand::MergeGroupsCommand(string option) {
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
			outputTypes["group"] = tempOutNames;
            outputTypes["count"] = tempOutNames;
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";	}
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("design");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["design"] = inputDir + it->second;		}
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
                
                it = parameters.find("count");
                //user has given a template file
                if(it != parameters.end()){
                    path = m->hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["count"] = inputDir + it->second;		}
                }
                
                it = parameters.find("fasta");
                //user has given a template file
                if(it != parameters.end()){
                    path = m->hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["fasta"] = inputDir + it->second;		}
                }
				
			}
			
			//check for required parameters
			designfile = validParameter.validFile(parameters, "design", true);
			if (designfile == "not open") { abort = true; }
			else if (designfile == "not found") {  				
				//if there is a current shared file, use it
				designfile = m->getDesignFile(); 
				if (designfile != "") { m->mothurOut("Using " + designfile + " as input file for the design parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current designfile and the design parameter is required."); m->mothurOutEndLine(); abort = true; }
			}else { m->setDesignFile(designfile); }	
			
			sharedfile = validParameter.validFile(parameters, "shared", true);
			if (sharedfile == "not open") { abort = true; sharedfile = ""; }
			else if (sharedfile == "not found") {  sharedfile = ""; }
			else { m->setSharedFile(sharedfile); }	
			
			groupfile = validParameter.validFile(parameters, "group", true);
			if (groupfile == "not open") { abort = true; groupfile = ""; }
			else if (groupfile == "not found") {  groupfile = ""; }
			else { m->setGroupFile(groupfile); }
            
            countfile = validParameter.validFile(parameters, "count", true);
            if (countfile == "not open") { abort = true; countfile = ""; }
            else if (countfile == "not found") {  countfile = ""; }
            else { m->setCountTableFile(countfile); }
            
            fastafile = validParameter.validFile(parameters, "fasta", true);
            if (fastafile == "not open") { abort = true; countfile = ""; }
            else if (fastafile == "not found") {  fastafile = ""; }
            else { m->setFastaFile(fastafile); }
            
			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  m->splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
			
			groups = validParameter.validFile(parameters, "groups", false);			
			if (groups == "not found") { groups = "all";  }
			m->splitAtDash(groups, Groups);
			m->setGroups(Groups);
            
            method = validParameter.validFile(parameters, "method", false);		if(method == "not found"){	method = "sum"; }
            
            if ((method != "sum") && (method != "average") && (method != "median")) { m->mothurOut(method + " is not a valid method. Options are sum, average and median. I will use sum."); m->mothurOutEndLine(); method = "sum"; }
            
            
            if ((groupfile != "") && (countfile != "")) {
                m->mothurOut("[ERROR]: you may only use one of the following: group or count."); m->mothurOutEndLine(); abort=true;
            }
			
			if ((sharedfile == "") && (groupfile == "") && (countfile == "")) {
				//give priority to group, then shared
				groupfile = m->getGroupFile(); 
				if (groupfile != "") {  m->mothurOut("Using " + groupfile + " as input file for the group parameter."); m->mothurOutEndLine(); }
				else { 
					sharedfile = m->getSharedFile(); 
					if (sharedfile != "") { m->mothurOut("Using " + sharedfile + " as input file for the shared parameter."); m->mothurOutEndLine(); }
					else { 
                        countfile = m->getCountTableFile();
                        if (countfile != "") { m->mothurOut("Using " + countfile + " as input file for the count parameter."); m->mothurOutEndLine(); }
                        else {
                            m->mothurOut("You have no current groupfile, countfile or sharedfile and one is required."); m->mothurOutEndLine(); abort = true;
                        }
					}
				}
			}
            
            if ((countfile == "") && (fastafile != "")) { m->mothurOut("[ERROR]: You may only use the fasta file with the count file, quitting."); m->mothurOutEndLine(); abort=true; }
            else if ((countfile != "") && (method == "average")) { m->mothurOut("You may not use the average method with the count file. I will use the sum method."); m->mothurOutEndLine(); method = "sum"; }
            else if ((countfile != "") && (method == "median") && (fastafile == "")) {
                fastafile = m->getFastaFile();
                if (fastafile != "") { m->mothurOut("Using " + fastafile + " as input file for the fasta parameter."); m->mothurOutEndLine(); }
                else {
                    m->mothurOut("[ERROR]: Fasta file is required with the median method and a count file so that sequences removed from your count table can also be removed from your fasta file to avoid downstream file mismatches, quitting.\n"); abort=true;
                }
            }
        
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "MergeGroupsCommand", "MergeGroupsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int MergeGroupsCommand::execute(){
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
	
		designMap = new DesignMap(designfile);
        
        if (method != "sum") {
            string defaultClass = designMap->getDefaultClass();
            vector<string> treatments = designMap->getCategory(defaultClass);
            set<int> numGroupsPerTreatment;
            for (int i = 0; i < treatments.size(); i++) {
                if (m->control_pressed) { break; }
                map<string, vector<string> > checkTreatments;
                vector<string> temp; temp.push_back(treatments[i]);
                checkTreatments[defaultClass] = temp;
                numGroupsPerTreatment.insert(designMap->getNumUnique(checkTreatments));
            }
            if (numGroupsPerTreatment.size() > 1) { m->mothurOut("[ERROR]: The median and average methods require you to have the same number of sequences in each treatment, quitting.\n"); delete designMap; return 0; }
        }

		if (groupfile != "") { processGroupFile(designMap); }
		if (sharedfile != "") { processSharedFile(designMap); }
        if (countfile != "") { processCountFile(designMap); }

		//reset groups parameter
		m->clearGroups();  
		delete designMap;
		
		if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); } return 0;}
		
		
		//set shared file as new current sharedfile
		string current = "";
		itTypes = outputTypes.find("shared");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setSharedFile(current); }
		}
		
		itTypes = outputTypes.find("group");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setGroupFile(current); }
		}
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "MergeGroupsCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************

int MergeGroupsCommand::process(vector<SharedRAbundVector*>& thisLookUp, ofstream& out){
	try {
        if (method == "average") {
            //create sharedRabundFloatVectors
            vector<SharedRAbundFloatVector*> temp = thisLookUp[0]->getSharedRAbundFloatVectors(thisLookUp);
            
            //follow code below
            map<string, SharedRAbundFloatVector> merged;
            map<string, SharedRAbundFloatVector>::iterator it;
            map<string, vector<int> > clearGroupAbunds;
            map<string, vector<int> >::iterator itAbunds;
            
            for (int i = 0; i < temp.size(); i++) {
                if (m->control_pressed) { return 0; }
                //what grouping does this group belong to
                string grouping = designMap->get(temp[i]->getGroup());
                if (grouping == "not found") { m->mothurOut("[ERROR]: " + temp[i]->getGroup() + " is not in your design file. Ignoring!"); m->mothurOutEndLine(); grouping = "NOTFOUND"; }
                else {
                    //do we already have a member of this grouping?
                    it = merged.find(grouping);
                    
                    if (it == merged.end()) { //nope, so create it
                        merged[grouping] = *temp[i];
                        merged[grouping].setGroup(grouping);
                        vector<int> temp;
                        clearGroupAbunds[grouping] = temp;
                    }
                }
            }
            
            for (int j = 0; j < temp[0]->getNumBins(); j++) {
                if (m->control_pressed) { return 0; }
                
                map<string, vector<int> > otusGroupAbunds = clearGroupAbunds;
                for (int i = 0; i < temp.size(); i++) {
                    
                    string grouping = designMap->get(temp[i]->getGroup());
                    if (grouping == "not found") { m->mothurOut("[ERROR]: " + temp[i]->getGroup() + " is not in your design file. Ignoring!"); m->mothurOutEndLine(); grouping = "NOTFOUND"; }
                    else {
                        otusGroupAbunds[grouping].push_back(temp[i]->getAbundance(j));
                    }
                }
                
                for (itAbunds = otusGroupAbunds.begin(); itAbunds != otusGroupAbunds.end(); itAbunds++) {
                    int abund = mergeAbund(itAbunds->second);
                    merged[itAbunds->first].set(j, abund, itAbunds->first);
                }
            }
            
            if (method == "median") {
                vector<SharedRAbundFloatVector*> temp2;
                for (it = merged.begin(); it != merged.end(); it++) {  temp2.push_back(&(it->second)); }
                eliminateZeroOTUS(temp2);
            }
            
            //print new file
            for (it = merged.begin(); it != merged.end(); it++) {
                if (!m->printedSharedHeaders) { (it->second).printHeaders(out); }
                out << (it->second).getLabel() << '\t' << it->first << '\t';
                (it->second).print(out);
            }
        }else {
            map<string, SharedRAbundVector> merged;
            map<string, SharedRAbundVector>::iterator it;
            map<string, vector<int> > clearGroupAbunds;
            map<string, vector<int> >::iterator itAbunds;
            
            for (int i = 0; i < thisLookUp.size(); i++) {
                if (m->control_pressed) { return 0; }
                //what grouping does this group belong to
                string grouping = designMap->get(thisLookUp[i]->getGroup());
                if (grouping == "not found") { m->mothurOut("[ERROR]: " + thisLookUp[i]->getGroup() + " is not in your design file. Ignoring!"); m->mothurOutEndLine(); grouping = "NOTFOUND"; }
                else {
                    //do we already have a member of this grouping?
                    it = merged.find(grouping);
                    
                    if (it == merged.end()) { //nope, so create it
                        merged[grouping] = *thisLookUp[i];
                        merged[grouping].setGroup(grouping);
                        vector<int> temp;
                        clearGroupAbunds[grouping] = temp;
                    }
                }
            }
            
            for (int j = 0; j < thisLookUp[0]->getNumBins(); j++) {
                if (m->control_pressed) { return 0; }
                
                map<string, vector<int> > otusGroupAbunds = clearGroupAbunds;
                for (int i = 0; i < thisLookUp.size(); i++) {
                    
                    string grouping = designMap->get(thisLookUp[i]->getGroup());
                    if (grouping == "not found") { m->mothurOut("[ERROR]: " + thisLookUp[i]->getGroup() + " is not in your design file. Ignoring!"); m->mothurOutEndLine(); grouping = "NOTFOUND"; }
                    else {
                        otusGroupAbunds[grouping].push_back(thisLookUp[i]->getAbundance(j));
                    }
                }
                
                for (itAbunds = otusGroupAbunds.begin(); itAbunds != otusGroupAbunds.end(); itAbunds++) {
                    int abund = mergeAbund(itAbunds->second);
                    merged[itAbunds->first].set(j, abund, itAbunds->first);
                }
            }
            
            if (method == "median") {
                vector<SharedRAbundVector*> temp;
                for (it = merged.begin(); it != merged.end(); it++) {  temp.push_back(&(it->second)); }
                eliminateZeroOTUS(temp);
            }
            
            //print new file
            for (it = merged.begin(); it != merged.end(); it++) {
                if (!m->printedSharedHeaders) { (it->second).printHeaders(out); }
                out << (it->second).getLabel() << '\t' << it->first << '\t';
                (it->second).print(out);
            }
        }
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "MergeGroupsCommand", "process");
		exit(1);
	}
}
//**********************************************************************************************************************

int MergeGroupsCommand::processSharedFile(DesignMap*& designMap){
	try {
		
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += m->hasPath(sharedfile);  }
        map<string, string> variables; 
        variables["[filename]"] = thisOutputDir + m->getRootName(m->getSimpleName(sharedfile));
        variables["[extension]"] = m->getExtension(sharedfile);
		string outputFileName = getOutputFileName("shared", variables);
        outputTypes["shared"].push_back(outputFileName); outputNames.push_back(outputFileName);
		
		ofstream out;
		m->openOutputFile(outputFileName, out);
		
		InputData input(sharedfile, "sharedfile");
		lookup = input.getSharedRAbundVectors();
		string lastLabel = lookup[0]->getLabel();
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;
		
		//as long as you are not at the end of the file or done wih the lines you want
		while((lookup[0] != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			
			if (m->control_pressed) {  out.close(); for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  } m->clearGroups();  delete designMap;  for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); } return 0; }
			
			if(allLines == 1 || labels.count(lookup[0]->getLabel()) == 1){			
				
				m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
				process(lookup, out);
				
				processedLabels.insert(lookup[0]->getLabel());
				userLabels.erase(lookup[0]->getLabel());
			}
			
			if ((m->anyLabelsToProcess(lookup[0]->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
				string saveLabel = lookup[0]->getLabel();
				
				for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }  
				lookup = input.getSharedRAbundVectors(lastLabel);
				m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
				
				process(lookup, out);
				
				processedLabels.insert(lookup[0]->getLabel());
				userLabels.erase(lookup[0]->getLabel());
				
				//restore real lastlabel to save below
				lookup[0]->setLabel(saveLabel);
			}
			
			lastLabel = lookup[0]->getLabel();
			//prevent memory leak
			for (int i = 0; i < lookup.size(); i++) {  delete lookup[i]; lookup[i] = NULL; }
			
			if (m->control_pressed) {  out.close(); m->clearGroups();   delete designMap;  for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); } return 0; }
			
			//get next line to process
			lookup = input.getSharedRAbundVectors();				
		}
		
		if (m->control_pressed) { out.close(); m->clearGroups();  delete designMap;  for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); }  return 0; }
		
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
			lookup = input.getSharedRAbundVectors(lastLabel);
			
			m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
			process(lookup, out);
			
			for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }
		}
		
		out.close();
		
				
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "MergeGroupsCommand", "processSharedFile");
		exit(1);
	}
}
//**********************************************************************************************************************

int MergeGroupsCommand::processGroupFile(DesignMap*& designMap){
	try {
		
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += m->hasPath(groupfile);  }
        map<string, string> variables; 
        variables["[filename]"] = thisOutputDir + m->getRootName(m->getSimpleName(groupfile));
        variables["[extension]"] = m->getExtension(groupfile);
		string outputFileName = getOutputFileName("group", variables);
		outputTypes["group"].push_back(outputFileName); outputNames.push_back(outputFileName);
		
		ofstream out;
		m->openOutputFile(outputFileName, out);
		
		//read groupfile
		GroupMap groupMap(groupfile);
		groupMap.readMap();
		
		//fill Groups - checks for "all" and for any typo groups
		SharedUtil* util = new SharedUtil();
		vector<string> nameGroups = groupMap.getNamesOfGroups();
		util->setGroups(Groups, nameGroups);
		delete util;
		
		vector<string> namesOfSeqs = groupMap.getNamesSeqs();
		bool error = false;
		
		for (int i = 0; i < namesOfSeqs.size(); i++) {
			
			if (m->control_pressed) { break; }
			
			string thisGroup = groupMap.getGroup(namesOfSeqs[i]);
			
			//are you in a group the user wants
			if (m->inUsersGroups(thisGroup, Groups)) {
				string thisGrouping = designMap->get(thisGroup);
				
				if (thisGrouping == "not found") { m->mothurOut("[ERROR]: " + namesOfSeqs[i] + " is from group " + thisGroup + " which is not in your design file, please correct."); m->mothurOutEndLine();  error = true; }
				else {
					out << namesOfSeqs[i] << '\t' << thisGrouping << endl;
				}
			}
		}
		
		if (error) { m->control_pressed = true; }

		out.close();
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "MergeGroupsCommand", "processGroupFile");
		exit(1);
	}
}
//**********************************************************************************************************************

int MergeGroupsCommand::processCountFile(DesignMap*& designMap){
    try {
        CountTable countTable;
        if (!countTable.testGroups(countfile)) { m->mothurOut("[ERROR]: your countfile contains no group information, please correct.\n"); m->control_pressed = true; return 0; }
        
        //read countTable
        countTable.readTable(countfile, true, false);
        
        //fill Groups - checks for "all" and for any typo groups
        SharedUtil util;
        vector<string> nameGroups = countTable.getNamesOfGroups();
        util.setGroups(Groups, nameGroups);
        
        vector<string> dnamesGroups = designMap->getNamesGroups();
        
        //sanity check
        bool error = false;
        if (nameGroups.size() == dnamesGroups.size()) { //at least there are the same number
            //is every group in counttable also in designmap
            for (int i = 0; i < nameGroups.size(); i++) {
                if (m->control_pressed) { break; }
                if (!m->inUsersGroups(nameGroups[i], dnamesGroups)) { error = true; break; }
            }
            
        }
        if (error) { m->mothurOut("[ERROR]: Your countfile does not contain the same groups as your design file, please correct\n"); m->control_pressed = true; return 0; }
        
        //user selected groups - remove some groups from table
        if (Groups.size() != nameGroups.size()) {
            for (int i = 0; i < nameGroups.size(); i++) {
                if (!m->inUsersGroups(nameGroups[i], Groups)) { countTable.removeGroup(nameGroups[i]); }
            }
        }
        //ask again in case order changed
        nameGroups = countTable.getNamesOfGroups();
        int numGroups = nameGroups.size();
        
        //create new table
        CountTable newTable;
        vector<string> treatments = designMap->getCategory();
        map<string, vector<int> > clearedMap;
        for (int i = 0; i < treatments.size(); i++) {
            newTable.addGroup(treatments[i]);
            vector<int> temp;
            clearedMap[treatments[i]] = temp;
        }
        treatments = newTable.getNamesOfGroups();
        
        set<string> namesToRemove;
        vector<string> namesOfSeqs = countTable.getNamesOfSeqs();
        for (int i = 0; i < namesOfSeqs.size(); i++) {
            
            if (m->control_pressed) { break; }
            
            vector<int> thisSeqsCounts = countTable.getGroupCounts(namesOfSeqs[i]);
            map<string, vector<int> > thisSeqsMap = clearedMap;
            
            for (int j = 0; j < numGroups; j++) {
                thisSeqsMap[designMap->get(nameGroups[j])].push_back(thisSeqsCounts[j]);
            }
        
            //create new counts for seq for new table
            vector<int> newCounts; int totalAbund = 0;
            for (int j = 0; j < treatments.size(); j++){
                int abund = mergeAbund(thisSeqsMap[treatments[j]]);
                newCounts.push_back(abund);  //order matters, add in count for each treatment in new table.
                totalAbund += abund;
            }
            
            //add seq to new table
            if(totalAbund == 0) {
                namesToRemove.insert(namesOfSeqs[i]);
            }else { newTable.push_back(namesOfSeqs[i], newCounts); }
        }
        
        if (error) { m->control_pressed = true; return 0; }
        
        //remove sequences zeroed out by median method
        if (namesToRemove.size() != 0) {
            //print names
            ofstream out;
            string accnosFile = "accnosFile.temp";
            m->openOutputFile(accnosFile, out);
            
            //output to .accnos file
            for (set<string>::iterator it = namesToRemove.begin(); it != namesToRemove.end(); it++) {
                if (m->control_pressed) {  out.close(); m->mothurRemove(accnosFile); return 0; }
                out << *it << endl;
            }
            out.close();

            //run remove.seqs
            string inputString = "accnos=" + accnosFile + ", fasta=" + fastafile;
            
            m->mothurOut("/******************************************/"); m->mothurOutEndLine();
            m->mothurOut("Running command: remove.seqs(" + inputString + ")"); m->mothurOutEndLine();
            m->mothurCalling = true;
            
            Command* removeCommand = new RemoveSeqsCommand(inputString);
            removeCommand->execute();
            
            map<string, vector<string> > filenames = removeCommand->getOutputFiles();
            
            delete removeCommand;
            m->mothurCalling = false;
            m->mothurOut("/******************************************/"); m->mothurOutEndLine();
            
            m->mothurRemove(accnosFile);
        }
    
        string thisOutputDir = outputDir;
        if (outputDir == "") {  thisOutputDir += m->hasPath(countfile);  }
        map<string, string> variables;
        variables["[filename]"] = thisOutputDir + m->getRootName(m->getSimpleName(countfile));
        variables["[extension]"] = m->getExtension(countfile);
        string outputFileName = getOutputFileName("count", variables);
        outputTypes["count"].push_back(outputFileName); outputNames.push_back(outputFileName);
        
        newTable.printTable(outputFileName);
        
        return 0;
        
    }
    catch(exception& e) {
        m->errorOut(e, "MergeGroupsCommand", "processCountFile");
        exit(1);
    }
}
//**********************************************************************************************************************

int MergeGroupsCommand::mergeAbund(vector<int> values){
    try {
        int abund = 0;
        
        if (method == "sum") {
            abund = m->sum(values);
        }else if (method == "average") {
            abund = m->average(values);
        }else if (method == "median") {
            abund = m->median(values);
        }else {
            m->mothurOut("[ERROR]: Invalid method. \n"); m->control_pressed = true; return 0;
        }
        
        return abund;
    }
    catch(exception& e) {
        m->errorOut(e, "MergeGroupsCommand", "mergeAbund");
        exit(1);
    }
}
//**********************************************************************************************************************
int MergeGroupsCommand::eliminateZeroOTUS(vector<SharedRAbundVector*>& thislookup) {
    try {
        
        vector<SharedRAbundVector*> newLookup;
        for (int i = 0; i < thislookup.size(); i++) {
            SharedRAbundVector* temp = new SharedRAbundVector();
            temp->setLabel(thislookup[i]->getLabel());
            temp->setGroup(thislookup[i]->getGroup());
            newLookup.push_back(temp);
        }
        
        //for each bin
        vector<string> newBinLabels;
        string snumBins = toString(thislookup[0]->getNumBins());
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
                //if there is a bin label use it otherwise make one
                string binLabel = "Otu";
                string sbinNumber = toString(i+1);
                if (sbinNumber.length() < snumBins.length()) {
                    int diff = snumBins.length() - sbinNumber.length();
                    for (int h = 0; h < diff; h++) { binLabel += "0"; }
                }
                binLabel += sbinNumber;
                if (i < m->currentSharedBinLabels.size()) {  binLabel = m->currentSharedBinLabels[i]; }
                
                newBinLabels.push_back(binLabel);
            }
        }
        
        for (int j = 0; j < thislookup.size(); j++) {  delete thislookup[j];  }
        thislookup.clear();
        
        thislookup = newLookup;
        m->currentSharedBinLabels = newBinLabels;
        
        return 0;
        
    }
    catch(exception& e) {
        m->errorOut(e, "MergeGroupsCommand", "eliminateZeroOTUS");
        exit(1);
    }
}
//**********************************************************************************************************************
int MergeGroupsCommand::eliminateZeroOTUS(vector<SharedRAbundFloatVector*>& thislookup) {
    try {
        
        vector<SharedRAbundFloatVector*> newLookup;
        for (int i = 0; i < thislookup.size(); i++) {
            SharedRAbundFloatVector* temp = new SharedRAbundFloatVector();
            temp->setLabel(thislookup[i]->getLabel());
            temp->setGroup(thislookup[i]->getGroup());
            newLookup.push_back(temp);
        }
        
        //for each bin
        vector<string> newBinLabels;
        string snumBins = toString(thislookup[0]->getNumBins());
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
                //if there is a bin label use it otherwise make one
                string binLabel = "Otu";
                string sbinNumber = toString(i+1);
                if (sbinNumber.length() < snumBins.length()) {
                    int diff = snumBins.length() - sbinNumber.length();
                    for (int h = 0; h < diff; h++) { binLabel += "0"; }
                }
                binLabel += sbinNumber;
                if (i < m->currentSharedBinLabels.size()) {  binLabel = m->currentSharedBinLabels[i]; }
                
                newBinLabels.push_back(binLabel);
            }
        }
        
        for (int j = 0; j < thislookup.size(); j++) {  delete thislookup[j];  }
        
        thislookup = newLookup;
        m->currentSharedBinLabels = newBinLabels;
        
        return 0;
        
    }
    catch(exception& e) {
        m->errorOut(e, "MergeGroupsCommand", "eliminateZeroOTUS");
        exit(1);
    }
}
//**********************************************************************************************************************



