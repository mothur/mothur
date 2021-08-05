//
//  listotucommand.cpp
//  Mothur
//
//  Created by Sarah Westcott on 5/15/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#include "listotuscommand.h"
#include "inputdata.h"

//**********************************************************************************************************************
vector<string> ListOtusCommand::setParameters(){	
	try {
        CommandParameter pshared("shared", "InputTypes", "", "", "SharedRel", "SharedRel", "none","accnos",false,false,true); parameters.push_back(pshared);
		CommandParameter prelabund("relabund", "InputTypes", "", "", "SharedRel", "SharedRel", "none","accnos",false,false); parameters.push_back(prelabund);
        CommandParameter plist("list", "InputTypes", "", "", "SharedRel", "SharedRel", "none","accnos",false,false); parameters.push_back(plist);
        CommandParameter pconstaxonomy("constaxonomy", "InputTypes", "", "", "SharedRel", "SharedRel", "none","accnos",false,false); parameters.push_back(pconstaxonomy);
        CommandParameter pgroups("groups", "String", "", "", "", "", "","",false,false); parameters.push_back(pgroups);
		CommandParameter plabel("label", "String", "", "", "", "", "","",false,false); parameters.push_back(plabel);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
        
        abort = false; calledHelp = false; allLines = true;
        
        vector<string> tempOutNames;
        outputTypes["accnos"] = tempOutNames;
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ListOtusCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string ListOtusCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The list.otus lists otu labels from shared, relabund, list or constaxonomy file. The results can be used by the get.otus to select specific otus with the output from classify.otu, otu.association, or corr.axes.\n";
		helpString += "The list.otulabels parameters are: shared, relabund, label and groups.\n";
		helpString += "The label parameter is used to analyze specific labels in your input.\n";
		helpString += "The groups parameter allows you to specify which of the groups you would like analyzed.\n";
		helpString += "The list.otulabels commmand should be in the following format: \n";
		helpString += "list.otulabels(shared=yourSharedFile, groups=yourGroup1-yourGroup2)\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "ListOtusCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string ListOtusCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "accnos") {  pattern = "[filename],[distance],accnos-[filename],accnos"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "ListOtusCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
ListOtusCommand::ListOtusCommand(string option) : Command()  {
	try {
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }
		
		else {
			OptionParser parser(option, setParameters());
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
            sharedfile = validParameter.validFile(parameters, "shared");
			if (sharedfile == "not open") { abort = true; }
			else if (sharedfile == "not found") { sharedfile = ""; }
			else { inputFileName = sharedfile; format = "sharedfile"; current->setSharedFile(sharedfile); }
			
			relabundfile = validParameter.validFile(parameters, "relabund");
			if (relabundfile == "not open") { abort = true; }
			else if (relabundfile == "not found") { relabundfile = ""; }
			else { inputFileName = relabundfile; format = "relabund"; current->setRelAbundFile(relabundfile); }
            
            listfile = validParameter.validFile(parameters, "list");
			if (listfile == "not open") { abort = true; }
			else if (listfile == "not found") { listfile = ""; }
			else { inputFileName = listfile; format = "list"; current->setListFile(listfile); }
            
            constaxonomy = validParameter.validFile(parameters, "constaxonomy");
            if (constaxonomy == "not open") { abort = true; }
            else if (constaxonomy == "not found") { constaxonomy = ""; }
            else { inputFileName = constaxonomy; format = "constaxonomy"; current->setConsTaxonomyFile(constaxonomy); }

            
            if ((relabundfile == "") && (sharedfile == "") && (constaxonomy == "") && (listfile== "")) {
				//is there are current file available for either of these?
				//give priority to shared, then relabund
				sharedfile = current->getSharedFile(); 
				if (sharedfile != "") {  inputFileName = sharedfile; format="sharedfile"; m->mothurOut("Using " + sharedfile + " as input file for the shared parameter.\n");  }
				else { 
					relabundfile = current->getRelAbundFile(); 
					if (relabundfile != "") {  inputFileName = relabundfile; format="relabund"; m->mothurOut("Using " + relabundfile + " as input file for the relabund parameter.\n");  }
					else { 
                        listfile = current->getListFile();
						if (listfile != "") {  inputFileName = listfile; format="list"; m->mothurOut("Using " + listfile + " as input file for the list parameter.\n");  }
                        else { 
                            constaxonomy = current->getConsTaxonomyFile();
                            if (constaxonomy != "") {  inputFileName = constaxonomy; format="constaxonomy"; m->mothurOut("Using " + constaxonomy + " as input file for the constaxonomy parameter.\n");  }
                            else {
                                m->mothurOut("No valid current files. You must provide a shared, list, relabund or constaxonomy file.\n");
                                abort = true;
                            }
                        }
					}
				}
			}
             
			if (outputdir == ""){ outputdir = util.hasPath(inputFileName);  }
            
            string groups = validParameter.valid(parameters, "groups");			
			if (groups == "not found") { groups = ""; }
			else { util.splitAtDash(groups, Groups); if (Groups.size() != 0) { if (Groups[0]== "all") { Groups.clear(); } } }
            
            string label = validParameter.valid(parameters, "label");			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  util.splitAtDash(label, labels);  allLines = false;  }
				else { allLines = true;  }
			}			
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "ListOtusCommand", "ListOtusCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int ListOtusCommand::execute(){
	try {
		
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
        
        InputData input(inputFileName, format, Groups);
        set<string> processedLabels;
        set<string> userLabels = labels;
        string lastLabel = "";
        
        if (format == "relabund") {
            SharedRAbundFloatVectors* lookup = util.getNextRelabund(input, allLines, userLabels, processedLabels, lastLabel);
            Groups = lookup->getNamesGroups();
            
            while (lookup != NULL) {
                
                if (m->getControl_pressed()) { delete lookup; break; }
                
                printList(lookup->getOTUNames(), lookup->getLabel()); delete lookup;
                
                lookup = util.getNextRelabund(input, allLines, userLabels, processedLabels, lastLabel);
            }

        }else if (format == "sharedfile") {
        
            SharedRAbundVectors* lookup = util.getNextShared(input, allLines, userLabels, processedLabels, lastLabel);
            Groups = lookup->getNamesGroups();
            
            while (lookup != NULL) {
                
                if (m->getControl_pressed()) { delete lookup; break; }
                
                printList(lookup->getOTUNames(), lookup->getLabel()); delete lookup;
                
                lookup = util.getNextShared(input, allLines, userLabels, processedLabels, lastLabel);
            }
            
        }else if (format == "list") {
            ListVector* list = util.getNextList(input, allLines, userLabels, processedLabels, lastLabel);
                   
            while (list != NULL) {
                       
                if (m->getControl_pressed()) { delete list; break; }
                       
                printList(list->getLabels(), list->getLabel()); delete list;
                      
                list = util.getNextList(input, allLines, userLabels, processedLabels, lastLabel);
            }
            
        }else if (format == "constaxonomy") { createList(constaxonomy); }
        
        if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) { util.mothurRemove(outputNames[i]); }  return 0; }
        
        //set relabund file as new current relabundfile
        string currentName = "";
        itTypes = outputTypes.find("accnos");
        if (itTypes != outputTypes.end()) {
            if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setAccnosFile(currentName); }
        }
        
        //output files created by command
		m->mothurOut("\nOutput File Names: \n"); 
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i] +"\n"); 	} m->mothurOutEndLine();
        return 0;
		
    }
	catch(exception& e) {
		m->errorOut(e, "ListOtusCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************

int ListOtusCommand::printList(vector<string> currentLabels, string distance){
	try {
        
        map<string, string> variables; 
        variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(inputFileName));
        variables["[distance]"] = distance;
        string outputFileName = getOutputFileName("accnos",variables);
        outputNames.push_back(outputFileName);  outputTypes["accnos"].push_back(outputFileName);
		ofstream out;
		util.openOutputFile(outputFileName, out);
        
        for (int i = 0; i < currentLabels.size(); i++) {  out << currentLabels[i] << endl;  }
        
        out.close();
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "ListOtusCommand", "printList");
		exit(1);
	}
}
//**********************************************************************************************************************

int ListOtusCommand::createList(string constaxFile){
    try {
        
        map<string, string> variables;
        variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(inputFileName));
        string outputFileName = getOutputFileName("accnos",variables);
        outputNames.push_back(outputFileName);  outputTypes["accnos"].push_back(outputFileName);
        ofstream out;
        util.openOutputFile(outputFileName, out);
        
        ifstream in;
        util.openInputFile(constaxFile, in);
        string otuLabel;
        
        //read headers
        string headers = util.getline(in);
        
        while(!in.eof()){
            
            if (m->getControl_pressed()) { break; }
            
            in >> otuLabel;
            string junk = util.getline(in); util.gobble(in);
            
            out << otuLabel << endl;
        }
        in.close();
        out.close();
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "ListOtusCommand", "createList");
        exit(1);
    }
}
//**********************************************************************************************************************

