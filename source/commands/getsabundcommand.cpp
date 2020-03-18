/*
 *  getsabundcommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 6/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "getsabundcommand.h"

//**********************************************************************************************************************
vector<string> GetSAbundCommand::setParameters(){	
	try {
		CommandParameter plist("list", "InputTypes", "", "", "LRSS", "LRSS", "none","sabund",false,false, true); parameters.push_back(plist);
        CommandParameter pcount("count", "InputTypes", "", "", "none", "none", "none","",false,false, false); parameters.push_back(pcount);
		CommandParameter prabund("rabund", "InputTypes", "", "", "LRSS", "LRSS", "none","sabund",false,false, true); parameters.push_back(prabund);		
		CommandParameter plabel("label", "String", "", "", "", "", "","",false,false); parameters.push_back(plabel);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "GetSAbundCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string GetSAbundCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The get.sabund command parameters is list, rabund, count and label.  list or rabund is required unless a valid current file exists.\n";
        helpString += "The count parameter allows you to provide a count file associated with your list file. If you clustered with a countfile the list file only contains the unique sequences and you will want to add the redundant counts into the sabund file, providing the count file allows you to do so.\n";
		helpString += "The label parameter allows you to select what distance levels you would like included in your .sabund file, and are separated by dashes.\n";
		helpString += "The get.sabund command should be in the following format: get.sabund(label=yourLabels).\n";
		helpString += "Example get.sabund().\n";
		helpString += "The default value for label is all labels in your inputfile.\n";
		helpString += "The get.sabund command outputs a .sabund file containing the labels you selected.\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "GetSAbundCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string GetSAbundCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "sabund")      {   pattern = "[filename],sabund";    }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "GetRAbundCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
GetSAbundCommand::GetSAbundCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["sabund"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "GetSAbundCommand", "GetSAbundCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
GetSAbundCommand::GetSAbundCommand(string option)  {
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
			map<string,string>::iterator it;
			
			ValidParameters validParameter;
			
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (!validParameter.isValidParameter(it->first, myArray, it->second)) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["sabund"] = tempOutNames;
			
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
				
				it = parameters.find("rabund");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["rabund"] = inputDir + it->second;		}
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
			if (listfile == "not open") { listfile = ""; abort = true; }
			else if (listfile == "not found") { listfile = ""; }
			else {  format = "list"; inputfile = listfile; current->setListFile(listfile); }
			
			rabundfile = validParameter.validFile(parameters, "rabund");
			if (rabundfile == "not open") { rabundfile = ""; abort = true; }	
			else if (rabundfile == "not found") { rabundfile = ""; }
			else {  format = "rabund"; inputfile = rabundfile; current->setRabundFile(rabundfile); }
			
            countfile = validParameter.validFile(parameters, "count");
			if (countfile == "not open") { countfile = ""; abort = true; }
			else if (countfile == "not found") { countfile = ""; }
			else {  current->setCountFile(countfile); }
            
						//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			label = validParameter.valid(parameters, "label");			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  util.splitAtDash(label, labels);  allLines = false;  }
				else { allLines = true;  }
			}
			
			if ((listfile == "") && (rabundfile == "")) { 
				//is there are current file available for any of these?
				//give priority to shared, then list, then rabund, then sabund
				//if there is a current shared file, use it
				listfile = current->getListFile(); 
				if (listfile != "") { inputfile = listfile; format = "list"; m->mothurOut("Using " + listfile + " as input file for the list parameter."); m->mothurOutEndLine(); }
				else { 
					rabundfile = current->getRabundFile(); 
					if (rabundfile != "") { inputfile = rabundfile; format = "rabund"; m->mothurOut("Using " + rabundfile + " as input file for the rabund parameter."); m->mothurOutEndLine(); }
					else { 
						m->mothurOut("No valid current files. You must provide a list or rabund file."); m->mothurOutEndLine(); 
						abort = true;
					}
				}
			}
			
            if ((countfile != "") && (listfile == "")) { m->mothurOut("[ERROR]: You can only use the count file with a list file, aborting.\n"); abort = true; }
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){	outputDir = util.hasPath(inputfile); 	}			
			
		}
	}
	catch(exception& e) {
		m->errorOut(e, "GetSAbundCommand", "GetSAbundCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int GetSAbundCommand::execute(){
	try {
		
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
        map<string, string> variables; 
		variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(inputfile));
		filename = getOutputFileName("sabund", variables);
		util.openOutputFile(filename, out);
		
        if (countfile != "") {
            processList(out);
        }else {
            InputData input(inputfile, format, nullVector);
            set<string> processedLabels;
            set<string> userLabels = labels;
            string lastLabel = "";
            
            SAbundVector* sabund = util.getNextSAbund(input, allLines, userLabels, processedLabels, lastLabel);
            
            while (sabund != NULL) {
                       
                if (m->getControl_pressed()) { delete sabund; break; }
                       
                sabund->print(out); delete sabund;
                      
                sabund = util.getNextSAbund(input, allLines, userLabels, processedLabels, lastLabel);
            }
		}
		out.close();
        
        if (m->getControl_pressed()) {  outputTypes.clear();  util.mothurRemove(filename);  return 0; }
		
		m->mothurOut("\nOutput File Names: \n"); 
		m->mothurOut(filename); m->mothurOutEndLine();	outputNames.push_back(filename); outputTypes["sabund"].push_back(filename);
		m->mothurOutEndLine();
		
		//set sabund file as new current sabundfile
		string currentName = "";
		itTypes = outputTypes.find("sabund");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setSabundFile(currentName); }
		}
		
		return 0;		
	}

	catch(exception& e) {
		m->errorOut(e, "GetSAbundCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
int GetSAbundCommand::processList(ofstream& out){
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
            SAbundVector sabund = rabund->getSAbundVector();
            sabund.print(out);
            delete rabund; delete list;
                  
            list = util.getNextList(input, allLines, userLabels, processedLabels, lastLabel);
        }
 
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "GetSAbundCommand", "processList");
		exit(1);
	}
}
//**********************************************************************************************************************

int GetSAbundCommand::createRabund(CountTable& ct, ListVector*& list, RAbundVector*& rabund){
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
		m->errorOut(e, "GetSAbundCommand", "createRabund");
		exit(1);
	}
    
}
//**********************************************************************************************************************




