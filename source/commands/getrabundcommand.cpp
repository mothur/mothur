/*
 *  getrabundcommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 6/2/09.
 *  Copyright 2009 Schloss Lab Umass Amherst. All rights reserved.
 *
 */

#include "getrabundcommand.h"

//**********************************************************************************************************************
vector<string> GetRAbundCommand::setParameters(){	
	try {
        CommandParameter pshared("shared", "InputTypes", "", "", "LRSS", "LRSS", "none","rabund",false,false, true); parameters.push_back(pshared);
		CommandParameter plist("list", "InputTypes", "", "", "LRSS", "LRSS", "none","rabund",false,false, true); parameters.push_back(plist);
        CommandParameter pcount("count", "InputTypes", "", "", "none", "none", "none","",false,false, false); parameters.push_back(pcount);
		CommandParameter psabund("sabund", "InputTypes", "", "", "LRSS", "LRSS", "none","rabund",false,false, true); parameters.push_back(psabund);		
		CommandParameter psorted("sorted", "Boolean", "", "T", "", "", "","",false,false); parameters.push_back(psorted);
		CommandParameter plabel("label", "String", "", "", "", "", "","",false,false); parameters.push_back(plabel);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pgroups("groups", "String", "", "", "", "", "","",false,false); parameters.push_back(pgroups);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
        
        abort = false; calledHelp = false; allLines = true;
        
        vector<string> tempOutNames;
        outputTypes["rabund"] = tempOutNames;
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "GetRAbundCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string GetRAbundCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The get.rabund command parameters are list, shared, sabund, count, label, groups and sorted.  shared, list or sabund parameters are required, unless you have valid current files.\n";
        helpString += "The count parameter allows you to provide a count file associated with your list file. If you clustered with a countfile the list file only contains the unique sequences and you will want to add the redundant counts into the rabund file, providing the count file allows you to do so.\n";
		helpString += "The label parameter allows you to select what distance levels you would like included in your .rabund file, and are separated by dashes.\n";
		helpString += "The sorted parameters allows you to print the rabund results sorted by abundance or not.  The default is sorted.\n";
		helpString += "The get.rabund command should be in the following format: get.rabund(label=yourLabels, sorted=yourSorted).\n";
		helpString += "Example get.rabund(sorted=F).\n";
		helpString += "The default value for label is all labels in your inputfile.\n";
		helpString += "The get.rabund command outputs a .rabund file containing the lines you selected.\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "GetRAbundCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string GetRAbundCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "rabund")      {   pattern = "[filename],rabund";    }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "GetRAbundCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
GetRAbundCommand::GetRAbundCommand(string option) : Command()  {
	try {
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }
		
		else {
			OptionParser parser(option, setParameters());
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			listfile = validParameter.validFile(parameters, "list");
			if (listfile == "not open") { listfile = ""; abort = true; }
			else if (listfile == "not found") { listfile = ""; }
			else {  format = "list"; inputfile = listfile; current->setListFile(listfile); }
			
            sharedfile = validParameter.validFile(parameters, "shared");
            if (sharedfile == "not open") { sharedfile = ""; abort = true; }
            else if (sharedfile == "not found") { sharedfile = ""; }
            else { format = "sharedfile"; inputfile = sharedfile; current->setSharedFile(sharedfile); }
            
			sabundfile = validParameter.validFile(parameters, "sabund");
			if (sabundfile == "not open") { sabundfile = ""; abort = true; }	
			else if (sabundfile == "not found") { sabundfile = ""; }
			else {  format = "sabund"; inputfile = sabundfile; current->setSabundFile(sabundfile); }
            
            countfile = validParameter.validFile(parameters, "count");
			if (countfile == "not open") { countfile = ""; abort = true; }
			else if (countfile == "not found") { countfile = ""; }
			else {  current->setCountFile(countfile); }
			
            string groups = validParameter.valid(parameters, "groups");
            if (groups == "not found") { groups = ""; }
            else {  util.splitAtDash(groups, Groups);  }

			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			string temp;
			temp = validParameter.valid(parameters, "sorted");			if (temp == "not found") { temp = "T"; }
			sorted = util.isTrue(temp);
			
			label = validParameter.valid(parameters, "label");			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  util.splitAtDash(label, labels);  allLines = false;  }
				else { allLines = true;  }
			}
			
			if ((listfile == "") && (sabundfile == "") && (sharedfile == "")) {
				//is there are current file available for any of these?
				//give priority to shared, then list, then rabund, then sabund
				//if there is a current shared file, use it
				listfile = current->getListFile(); 
				if (listfile != "") { inputfile = listfile; format = "list"; m->mothurOut("Using " + listfile + " as input file for the list parameter.\n");}
				else {
                    sharedfile = current->getSharedFile();
                    if (sharedfile != "") { inputfile = sharedfile; format = "sharedfile"; m->mothurOut("Using " + sharedfile + " as input file for the shared parameter.\n"); }
                    else {
                        sabundfile = current->getSabundFile();
                        if (sabundfile != "") { inputfile = sabundfile; format = "sabund"; m->mothurOut("Using " + sabundfile + " as input file for the sabund parameter.\n"); }
                        else { m->mothurOut("No valid current files. You must provide a shared, list or sabund file.\n"); abort = true; }
                    }
				}
			}
			
			if ((countfile != "") && (listfile == "")) { m->mothurOut("[ERROR]: You can only use the count file with a list file, aborting.\n"); abort = true; }
            
			 
					if (outputdir == ""){    outputdir = util.hasPath(inputfile); 	}			
		}
	}
	catch(exception& e) {
		m->errorOut(e, "GetRAbundCommand", "GetRAbundCommand");
		exit(1);
	}			
}
//**********************************************************************************************************************

int GetRAbundCommand::execute(){
	try {
	
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
        
        map<string, string> variables; 
        variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(inputfile));
		filename = getOutputFileName("rabund", variables);
		util.openOutputFile(filename, out);
        
        if (countfile != "") {
            processList(out);
        }else {
            InputData input(inputfile, format, nullVector);
            set<string> processedLabels;
            set<string> userLabels = labels;
            string lastLabel = "";
            
            RAbundVector* rabund = util.getNextRAbund(input, allLines, userLabels, processedLabels, lastLabel);
            
            while (rabund != NULL) {
                       
                if (m->getControl_pressed()) { delete rabund; break; }
                       
                if(sorted)      {   rabund->print(out);                 }
                else            {   rabund->nonSortedPrint(out);        }
                delete rabund;
                      
                rabund = util.getNextRAbund(input, allLines, userLabels, processedLabels, lastLabel);
            }
        }
        out.close();
        
        if (m->getControl_pressed()) {  outputTypes.clear();  util.mothurRemove(filename);  return 0; }
        
		m->mothurOut("\nOutput File Names:\n"+filename+"\n\n");
        outputNames.push_back(filename); outputTypes["rabund"].push_back(filename);
		
		//set rabund file as new current rabundfile
		string currentName = "";
		itTypes = outputTypes.find("rabund");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setRabundFile(currentName); }
		}
		
		return 0;		
	}
	catch(exception& e) {
		m->errorOut(e, "GetRAbundCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
int GetRAbundCommand::processList(ofstream& out){
	try {
        CountTable ct;
        ct.readTable(countfile, false, false);
        
        InputData input(inputfile, format, nullVector);
        set<string> processedLabels;
        set<string> userLabels = labels;
        string lastLabel = "";
        
        ListVector* list = util.getNextList(input, allLines, userLabels, processedLabels, lastLabel);
        
        if (m->getControl_pressed()) {  delete list;  return 0; }
        
        while (list != NULL) {
                   
            if (m->getControl_pressed()) { delete list; break; }
                   
            RAbundVector* rabund = new RAbundVector();
            createRabund(ct, list, rabund);
            
            if(sorted)  {   rabund->print(out);              }
            else        {    rabund->nonSortedPrint(out);    }
            
            delete rabund; delete list;
                  
            list = util.getNextList(input, allLines, userLabels, processedLabels, lastLabel);
        }
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "GetRAbundCommand", "processList");
		exit(1);
	}
}
//**********************************************************************************************************************

int GetRAbundCommand::createRabund(CountTable& ct, ListVector*& list, RAbundVector*& rabund){
    try {
        
        rabund->setLabel(list->getLabel());
        for(int i = 0; i < list->getNumBins(); i++) {
            if (m->getControl_pressed()) { return 0; }
            vector<string> binNames;
            string bin = list->get(i);
            util.splitAtComma(bin, binNames);
            int total = 0;
            for (int j = 0; j < binNames.size(); j++) {
                total += ct.getNumSeqs(binNames[j]);
            }
            rabund->push_back(total);
        }
        
        return 0;
    }
    catch(exception& e) {
		m->errorOut(e, "GetRAbundCommand", "createRabund");
		exit(1);
	}
}
//**********************************************************************************************************************


