/*
 *  getlistcountcommand.cpp
 *  Mothur
 *
 *  Created by westcott on 10/12/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */

#include "getlistcountcommand.h"

//**********************************************************************************************************************
vector<string> GetListCountCommand::setParameters(){	
	try {
		CommandParameter plist("list", "InputTypes", "", "", "none", "none", "none","otu",false,true, true); parameters.push_back(plist);
		CommandParameter plabel("label", "String", "", "", "", "", "","",false,false); parameters.push_back(plabel);
		CommandParameter parasort("sort", "Multiple", "name-otu", "otu", "", "", "","",false,false); parameters.push_back(parasort);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
        abort = false; calledHelp = false; allLines = true;
        
        vector<string> tempOutNames;
        outputTypes["otu"] = tempOutNames;
        
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "GetListCountCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string GetListCountCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The get.otulist command parameters are list, sort and label.  list is required, unless you have a valid current list file.\n";
		helpString += "The label parameter allows you to select what distance levels you would like a output files created for, and are separated by dashes.\n";
		helpString += "The sort parameter allows you to select how you want the output displayed. Options are otu and name.\n";
		helpString += "If otu is selected the output will be otu number followed by the list of names in that otu.\n";
		helpString += "If name is selected the output will be a sequence name followed by its otu number.\n";
		helpString += "The get.otulist command should be in the following format: get.otulist(list=yourlistFile, label=yourLabels).\n";
		helpString += "Example get.otulist(list=amazon.fn.list, label=0.10).\n";
		helpString += "The default value for label is all lines in your inputfile.\n";
		helpString += "The get.otulist command outputs a .otu file for each distance you specify listing the bin number and the names of the sequences in that bin.\n";
		
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "GetListCountCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string GetListCountCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "otu") {  pattern = "[filename],[tag],otu"; } 
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "GetListCountCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
GetListCountCommand::GetListCountCommand(string option)  {
	try {
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }
		
		else {
			OptionParser parser(option, setParameters());
			map<string, string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){	outputDir = "";		}
			
			//check for required parameters
			listfile = validParameter.validFile(parameters, "list");
			if (listfile == "not found")  { 				
				listfile = current->getListFile(); 
				if (listfile != "") { m->mothurOut("Using " + listfile + " as input file for the list parameter.\n"); }
				else { 	m->mothurOut("You have no current list file and the list parameter is required.\n");  abort = true; }
			}
			else if (listfile == "not open") { abort = true; }	
			else { current->setListFile(listfile); }
			
		
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			sort = validParameter.valid(parameters, "sort");	  if (sort == "not found") { sort = "otu"; }
			if ((sort != "otu") && (sort != "name")) { m->mothurOut( sort + " is not a valid sort option. Options are otu and name. I will use otu.\n");  sort = "otu"; }
			
			label = validParameter.valid(parameters, "label");			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  util.splitAtDash(label, labels);  allLines = false;  }
				else { allLines = true;  }
			}
		}
	}
	catch(exception& e) {
		m->errorOut(e, "GetListCountCommand", "GetListCountCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int GetListCountCommand::execute(){
	try {
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
        
        InputData input(listfile, "list", nullVector);
        set<string> processedLabels;
        set<string> userLabels = labels;
        string lastLabel = "";
        
        ListVector* list = util.getNextList(input, allLines, userLabels, processedLabels, lastLabel);
               
        while (list != NULL) {
                   
            if (m->getControl_pressed()) { delete list; break; }
                   
            process(list); delete list;
                  
            list = util.getNextList(input, allLines, userLabels, processedLabels, lastLabel);
        }
        
		if (m->getControl_pressed()) { return 0; }
        
		m->mothurOut("\nOutput File Names: \n"); 
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i] +"\n"); 	} m->mothurOutEndLine();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "GetListCountCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************
//return 1 if error, 0 otherwise
void GetListCountCommand::process(ListVector* list) {
	try {
		string binnames;
		if (outputDir == "") { outputDir += util.hasPath(listfile); }
        map<string, string> variables; 
		variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(listfile));
        variables["[tag]"] = list->getLabel();
		string outputFileName = getOutputFileName("otu", variables);
		
		util.openOutputFile(outputFileName, out);
		outputNames.push_back(outputFileName); outputTypes["otu"].push_back(outputFileName);
		
		m->mothurOut(list->getLabel()); m->mothurOutEndLine();
		
		//for each bin in the list vector
        vector<string> binLabels = list->getLabels();
		for (int i = 0; i < list->getNumBins(); i++) {
			if (m->getControl_pressed()) { break; }
			
			binnames = list->get(i);
			
			if (sort == "otu") {
				out << binLabels[i] << '\t' << binnames << endl;
			}else{ //sort = name
				vector<string> names;
				util.splitAtComma(binnames, names);
				
				for (int j = 0; j < names.size(); j++) {
					out << names[j] << '\t' << binLabels[i] << endl;
				}
			}
		}
		
		out.close();
	}
	catch(exception& e) {
		m->errorOut(e, "GetListCountCommand", "process");
		exit(1);
	}
}
//**********************************************************************************************************************


