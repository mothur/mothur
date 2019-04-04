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
GetListCountCommand::GetListCountCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["otu"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "GetListCountCommand", "GetListCountCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
GetListCountCommand::GetListCountCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		allLines = true;
				
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string, string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			map<string, string>::iterator it;
			
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (!validParameter.isValidParameter(it->first, myArray, it->second)) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["otu"] = tempOutNames;
			
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
			}

			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){	outputDir = "";		}
			
			//check for required parameters
			listfile = validParameter.validFile(parameters, "list");
			if (listfile == "not found")  { 				
				listfile = current->getListFile(); 
				if (listfile != "") { m->mothurOut("Using " + listfile + " as input file for the list parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current list file and the list parameter is required."); m->mothurOutEndLine(); abort = true; }
			}
			else if (listfile == "not open") { abort = true; }	
			else { current->setListFile(listfile); }
			
		
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			sort = validParameter.valid(parameters, "sort");	  if (sort == "not found") { sort = "otu"; }
			if ((sort != "otu") && (sort != "name")) { m->mothurOut( sort + " is not a valid sort option. Options are otu and name. I will use otu."); m->mothurOutEndLine(); sort = "otu"; }
			
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
		
		input = new InputData(listfile, "list", nullVector);
		list = input->getListVector();
		string lastLabel = list->getLabel();

		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;
		
		if (m->getControl_pressed()) { delete input; delete list; for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]);	} return 0;  }
		
		while((list != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			
			if(allLines == 1 || labels.count(list->getLabel()) == 1){
			
				process(list);
				
				if (m->getControl_pressed()) { delete input; delete list; for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]);	} return 0;  }
							
				processedLabels.insert(list->getLabel());
				userLabels.erase(list->getLabel());
			}
			
			if ((util.anyLabelsToProcess(list->getLabel(), userLabels, "") ) && (processedLabels.count(lastLabel) != 1)) {
				string saveLabel = list->getLabel();
				
				delete list;
				list = input->getListVector(lastLabel);
				
				process(list);
				
				if (m->getControl_pressed()) { delete input; delete list; for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]);	} return 0;  }

													
				processedLabels.insert(list->getLabel());
				userLabels.erase(list->getLabel());
				
				//restore real lastlabel to save below
				list->setLabel(saveLabel);
			}
			
			lastLabel = list->getLabel();			
			
			delete list;
			list = input->getListVector();
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
				m->mothurOut(". Please refer to " + lastLabel + ".");  m->mothurOutEndLine();
			}
		}
		
		//run last label if you need to
		if (needToRun )  {
			if (list != NULL) {		delete list;	}
			list = input->getListVector(lastLabel);
				
			process(list);	
			
			if (m->getControl_pressed()) { delete input; delete list; for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]);	} return 0;  }
			
			delete list;  
		}
		
		delete input;
		
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


