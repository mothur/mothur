//
//  listotucommand.cpp
//  Mothur
//
//  Created by Sarah Westcott on 5/15/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#include "listotulabelscommand.h"
#include "inputdata.h"

//**********************************************************************************************************************
vector<string> ListOtuLabelsCommand::setParameters(){	
	try {
        CommandParameter pshared("shared", "InputTypes", "", "", "SharedRel", "SharedRel", "none","otulabels",false,false,true); parameters.push_back(pshared);
		CommandParameter prelabund("relabund", "InputTypes", "", "", "SharedRel", "SharedRel", "none","otulabels",false,false); parameters.push_back(prelabund);
        CommandParameter plist("list", "InputTypes", "", "", "SharedRel", "SharedRel", "none","otulabels",false,false); parameters.push_back(plist);
        CommandParameter pgroups("groups", "String", "", "", "", "", "","",false,false); parameters.push_back(pgroups);
		CommandParameter plabel("label", "String", "", "", "", "", "","",false,false); parameters.push_back(plabel);
        //every command must have inputdir and outputdir.  This allows mothur users to redirect input and output files.
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ListOtuLabelsCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string ListOtuLabelsCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The list.otulabels lists otu labels from shared, relabund or list file. The results can be used by the get.otulabels to select specific otus with the output from classify.otu, otu.association, or corr.axes.\n";
		helpString += "The list.otulabels parameters are: shared, relabund, label and groups.\n";
		helpString += "The label parameter is used to analyze specific labels in your input.\n";
		helpString += "The groups parameter allows you to specify which of the groups you would like analyzed.\n";
		helpString += "The list.otulabels commmand should be in the following format: \n";
		helpString += "list.otulabels(shared=yourSharedFile, groups=yourGroup1-yourGroup2)\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "ListOtuLabelsCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string ListOtuLabelsCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "otulabels") {  pattern = "[filename],[distance],otulabels"; } 
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "ListOtuLabelsCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
ListOtuLabelsCommand::ListOtuLabelsCommand(){	
	try {
		abort = true; calledHelp = true;
		setParameters();
        vector<string> tempOutNames;
		outputTypes["otulabels"] = tempOutNames; 
	}
	catch(exception& e) {
		m->errorOut(e, "ListOtuLabelsCommand", "ListOtuLabelsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
ListOtuLabelsCommand::ListOtuLabelsCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		allLines = 1;
        
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			//valid paramters for this command
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			map<string,string>::iterator it;
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
                
                //edit file types below to include only the types you added as parameters
                
				string path;
                it = parameters.find("relabund");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["relabund"] = inputDir + it->second;		}
				}
                
                it = parameters.find("shared");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["shared"] = inputDir + it->second;		}
				}
                
                it = parameters.find("list");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["list"] = inputDir + it->second;		}
				}
            }
            
            vector<string> tempOutNames;
            outputTypes["otulabels"] = tempOutNames; 
            
 			//check for parameters
            sharedfile = validParameter.validFile(parameters, "shared", true);
			if (sharedfile == "not open") { abort = true; }
			else if (sharedfile == "not found") { sharedfile = ""; }
			else { inputFileName = sharedfile; format = "sharedfile"; m->setSharedFile(sharedfile); }
			
			relabundfile = validParameter.validFile(parameters, "relabund", true);
			if (relabundfile == "not open") { abort = true; }
			else if (relabundfile == "not found") { relabundfile = ""; }
			else { inputFileName = relabundfile; format = "relabund"; m->setRelAbundFile(relabundfile); }
            
            listfile = validParameter.validFile(parameters, "list", true);
			if (listfile == "not open") { abort = true; }
			else if (listfile == "not found") { listfile = ""; }
			else { inputFileName = listfile; format = "list"; m->setListFile(listfile); }

            
            if ((relabundfile == "") && (sharedfile == "") && (listfile== "")) { 
				//is there are current file available for either of these?
				//give priority to shared, then relabund
				sharedfile = m->getSharedFile(); 
				if (sharedfile != "") {  inputFileName = sharedfile; format="sharedfile"; m->mothurOut("Using " + sharedfile + " as input file for the shared parameter."); m->mothurOutEndLine(); }
				else { 
					relabundfile = m->getRelAbundFile(); 
					if (relabundfile != "") {  inputFileName = relabundfile; format="relabund"; m->mothurOut("Using " + relabundfile + " as input file for the relabund parameter."); m->mothurOutEndLine(); }
					else { 
                        listfile = m->getListFile();
						if (listfile != "") {  inputFileName = listfile; format="list"; m->mothurOut("Using " + listfile + " as input file for the list parameter."); m->mothurOutEndLine(); }
                        else { 
                            m->mothurOut("No valid current files. You must provide a shared, list or relabund."); m->mothurOutEndLine(); 
                            abort = true;
                        }
					}
				}
			}
            
            //if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = m->hasPath(inputFileName); //if user entered a file with a path then preserve it	
			}
            
            string groups = validParameter.validFile(parameters, "groups", false);			
			if (groups == "not found") { groups = ""; }
			else { m->splitAtDash(groups, Groups); }
			m->setGroups(Groups);
            
            string label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  m->splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}			
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "ListOtuLabelsCommand", "ListOtuLabelsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int ListOtuLabelsCommand::execute(){
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
        
        InputData input(inputFileName, format);
        
        if (format == "relabund") {
            vector<SharedRAbundFloatVector*> lookup = input.getSharedRAbundFloatVectors();
            string lastLabel = lookup[0]->getLabel();
            
            //if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
            set<string> processedLabels;
            set<string> userLabels = labels;
            
            //as long as you are not at the end of the file or done wih the lines you want
            while((lookup[0] != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
                
                if (m->control_pressed) { for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }  for (int i = 0; i < outputNames.size(); i++) { m->mothurRemove(outputNames[i]); } return 0; }
                
                if(allLines == 1 || labels.count(lookup[0]->getLabel()) == 1){			
                    
                    m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
                    
                    createList(lookup);
                    
                    processedLabels.insert(lookup[0]->getLabel());
                    userLabels.erase(lookup[0]->getLabel());
                }
                
                if ((m->anyLabelsToProcess(lookup[0]->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
                    string saveLabel = lookup[0]->getLabel();
                    
                    for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }  
                    lookup = input.getSharedRAbundFloatVectors(lastLabel);
                    m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
                    
                    createList(lookup);
                    
                    processedLabels.insert(lookup[0]->getLabel());
                    userLabels.erase(lookup[0]->getLabel());
                    
                    //restore real lastlabel to save below
                    lookup[0]->setLabel(saveLabel);
                }
                
                lastLabel = lookup[0]->getLabel();
                //prevent memory leak
                for (int i = 0; i < lookup.size(); i++) {  delete lookup[i]; lookup[i] = NULL; }
                
                if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) { m->mothurRemove(outputNames[i]); }  return 0; }
                
                //get next line to process
                lookup = input.getSharedRAbundFloatVectors();				
            }
            
            if (m->control_pressed) {  for (int i = 0; i < outputNames.size(); i++) { m->mothurRemove(outputNames[i]); }  return 0; }
            
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
                lookup = input.getSharedRAbundFloatVectors(lastLabel);
                
                m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
                
                createList(lookup);
                
                for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }
            }
        }else if (format == "sharedfile") {
            
            vector<SharedRAbundVector*> lookup = input.getSharedRAbundVectors();
            string lastLabel = lookup[0]->getLabel();
            
            //if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
            set<string> processedLabels;
            set<string> userLabels = labels;
            
            //as long as you are not at the end of the file or done wih the lines you want
            while((lookup[0] != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
                
                if (m->control_pressed) { for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }  for (int i = 0; i < outputNames.size(); i++) { m->mothurRemove(outputNames[i]); } return 0; }
                
                if(allLines == 1 || labels.count(lookup[0]->getLabel()) == 1){			
                    
                    m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
                    
                    createList(lookup);
                    
                    processedLabels.insert(lookup[0]->getLabel());
                    userLabels.erase(lookup[0]->getLabel());
                }
                
                if ((m->anyLabelsToProcess(lookup[0]->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
                    string saveLabel = lookup[0]->getLabel();
                    
                    for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }  
                    lookup = input.getSharedRAbundVectors(lastLabel);
                    m->mothurOut(lookup[0]->getLabel()); m->mothurOutEndLine();
                    
                    createList(lookup);
                    
                    processedLabels.insert(lookup[0]->getLabel());
                    userLabels.erase(lookup[0]->getLabel());
                    
                    //restore real lastlabel to save below
                    lookup[0]->setLabel(saveLabel);
                }
                
                lastLabel = lookup[0]->getLabel();
                //prevent memory leak
                for (int i = 0; i < lookup.size(); i++) {  delete lookup[i]; lookup[i] = NULL; }
                
                if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) { m->mothurRemove(outputNames[i]); }  return 0; }
                
                //get next line to process
                lookup = input.getSharedRAbundVectors();				
            }
            
            if (m->control_pressed) {  for (int i = 0; i < outputNames.size(); i++) { m->mothurRemove(outputNames[i]); }  return 0; }
            
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
                
                createList(lookup);
                
                for (int i = 0; i < lookup.size(); i++) {  delete lookup[i];  }
            }
        }else {
            ListVector* list = input.getListVector();
            string lastLabel = list->getLabel();
            
            //if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
            set<string> processedLabels;
            set<string> userLabels = labels;
            
            //as long as you are not at the end of the file or done wih the lines you want
            while((list != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
                
                if (m->control_pressed) { delete list;  for (int i = 0; i < outputNames.size(); i++) { m->mothurRemove(outputNames[i]); } return 0; }
                
                if(allLines == 1 || labels.count(list->getLabel()) == 1){			
                    
                    m->mothurOut(list->getLabel()); m->mothurOutEndLine();
                    
                    createList(list);
                    
                    processedLabels.insert(list->getLabel());
                    userLabels.erase(list->getLabel());
                }
                
                if ((m->anyLabelsToProcess(list->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
                    string saveLabel = list->getLabel();
                    
                    delete list; 
                    list = input.getListVector(lastLabel);
                    m->mothurOut(list->getLabel()); m->mothurOutEndLine();
                    
                    createList(list);
                    
                    processedLabels.insert(list->getLabel());
                    userLabels.erase(list->getLabel());
                    
                    //restore real lastlabel to save below
                    list->setLabel(saveLabel);
                }
                
                lastLabel = list->getLabel();
                //prevent memory leak
                delete list; list = NULL;
                
                if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) { m->mothurRemove(outputNames[i]); }  return 0; }
                
                //get next line to process
                list = input.getListVector();				
            }
            
            if (m->control_pressed) {  for (int i = 0; i < outputNames.size(); i++) { m->mothurRemove(outputNames[i]); }  return 0; }
            
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
                delete list;  
                list = input.getListVector(lastLabel);
                
                m->mothurOut(list->getLabel()); m->mothurOutEndLine();
                
                createList(list);
                
                delete list;
            }
        }
        
        if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) { m->mothurRemove(outputNames[i]); }  return 0; }
        
        //output files created by command
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
        return 0;
		
    }
	catch(exception& e) {
		m->errorOut(e, "ListOtuLabelsCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************

int ListOtuLabelsCommand::createList(vector<SharedRAbundVector*>& lookup){
	try {
        
        map<string, string> variables; 
        variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(inputFileName));
        variables["[distance]"] = lookup[0]->getLabel();
        string outputFileName = getOutputFileName("otulabels",variables);
        outputNames.push_back(outputFileName);  outputTypes["otulabels"].push_back(outputFileName);
		ofstream out;
		m->openOutputFile(outputFileName, out);
        
        for (int i = 0; i < m->currentSharedBinLabels.size(); i++) {  out << m->currentSharedBinLabels[i] << endl;  }
        
        out.close();
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "ListOtuLabelsCommand", "createList");
		exit(1);
	}
}

//**********************************************************************************************************************

int ListOtuLabelsCommand::createList(vector<SharedRAbundFloatVector*>& lookup){
	try {
        map<string, string> variables; 
        variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(inputFileName));
        variables["[distance]"] = lookup[0]->getLabel();
        string outputFileName = getOutputFileName("otulabels",variables);
        outputNames.push_back(outputFileName);  outputTypes["accnos"].push_back(outputFileName);
		ofstream out;
		m->openOutputFile(outputFileName, out);
        
        for (int i = 0; i < m->currentSharedBinLabels.size(); i++) {  out << m->currentSharedBinLabels[i] << endl;  }
        
        out.close();
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "ListOtuLabelsCommand", "createList");
		exit(1);
	}
}
//**********************************************************************************************************************
int ListOtuLabelsCommand::createList(ListVector*& list){
	try {
        map<string, string> variables; 
        variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(inputFileName));
        variables["[distance]"] = list->getLabel();
        string outputFileName = getOutputFileName("otulabels",variables);
        outputNames.push_back(outputFileName);  outputTypes["accnos"].push_back(outputFileName);
		ofstream out;
		m->openOutputFile(outputFileName, out);
        
        vector<string> binLabels = list->getLabels();
        for (int i = 0; i < binLabels.size(); i++) {  out << binLabels[i] << endl;  }

        out.close();
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "ListOtuLabelsCommand", "createList");
		exit(1);
	}
}

//**********************************************************************************************************************

