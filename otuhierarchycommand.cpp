/*
 *  otuhierarchycommand.cpp
 *  Mothur
 *
 *  Created by westcott on 1/19/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "otuhierarchycommand.h"

//**********************************************************************************************************************
vector<string> OtuHierarchyCommand::setParameters(){	
	try {
		CommandParameter poutput("output", "Multiple", "name-number", "name", "", "", "","",false,false); parameters.push_back(poutput);
		CommandParameter plist("list", "InputTypes", "", "", "none", "none", "none","otuheirarchy",false,true,true); parameters.push_back(plist);
		CommandParameter plabel("label", "String", "", "", "", "", "","",false,false); parameters.push_back(plabel);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "OtuHierarchyCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string OtuHierarchyCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The otu.hierarchy command is used to see how otus relate at two distances. \n";
		helpString += "The otu.hierarchy command parameters are list, label and output.  list and label parameters are required. \n";
		helpString += "The output parameter allows you to output the names of the sequence in the OTUs or the OTU numbers. Options are name and number, default is name. \n";
		helpString += "The otu.hierarchy command should be in the following format: \n";
		helpString += "otu.hierarchy(list=yourListFile, label=yourLabels).\n";
		helpString += "Example otu.hierarchy(list=amazon.fn.list, label=0.01-0.03).\n";
		helpString += "The otu.hierarchy command outputs a .otu.hierarchy file which is described on the wiki.\n";
		helpString += "Note: No spaces between parameter labels (i.e. list), '=' and parameters (i.e.yourListFile).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "OtuHierarchyCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string OtuHierarchyCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "otuheirarchy") {  pattern = "[filename],[distance1],[tag],[distance2],otu.hierarchy"; } 
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "OtuHierarchyCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
OtuHierarchyCommand::OtuHierarchyCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["otuheirarchy"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "OtuHierarchyCommand", "OtuHierarchyCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
OtuHierarchyCommand::OtuHierarchyCommand(string option) {
	try {
		abort = false; calledHelp = false;   
		
		//allow user to run help
		if(option == "help") {  help(); abort = true; calledHelp = true; }
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
			outputTypes["otuheirarchy"] = tempOutNames;
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("list");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["list"] = inputDir + it->second;		}
				}
			}

			listFile = validParameter.validFile(parameters, "list", true);
			if (listFile == "not found") { 
				listFile = m->getListFile(); 
				if (listFile != "") {  m->mothurOut("Using " + listFile + " as input file for the list parameter."); m->mothurOutEndLine(); }
				else { 
					m->mothurOut("No valid current list file. You must provide a list file."); m->mothurOutEndLine(); 
					abort = true;
				}
			}else if (listFile == "not open") { abort = true; }	
			else { m->setListFile(listFile); }
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += m->hasPath(listFile); //if user entered a file with a path then preserve it	
			}
			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { m->mothurOut("label is a required parameter for the otu.hierarchy command."); m->mothurOutEndLine(); abort = true; }
			else { 
				m->splitAtDash(label, labels);
				if (labels.size() != 2) { m->mothurOut("You must provide 2 labels."); m->mothurOutEndLine(); abort = true; }
			}	
			
			output = validParameter.validFile(parameters, "output", false);			if (output == "not found") { output = "name"; }
			
			if ((output != "name") && (output != "number")) { m->mothurOut("output options are name and number. I will use name."); m->mothurOutEndLine(); output = "name"; }
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "OtuHierarchyCommand", "OtuHierarchyCommand");
		exit(1);
	}			
}
//**********************************************************************************************************************

int OtuHierarchyCommand::execute(){
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		//get listvectors that correspond to labels requested, (or use smart distancing to get closest listvector)
		vector<ListVector> lists = getListVectors();
		
		if (m->control_pressed) { outputTypes.clear(); return 0; }
		
		//determine which is little and which is big, putting little first
		if (lists.size() == 2) {
			//if big is first swap them
			if (lists[0].getNumBins() < lists[1].getNumBins()) {
				reverse(lists.begin(), lists.end());
			}
		}else{
			m->mothurOut("error getting listvectors, unable to read 2 different vectors, check your label inputs."); m->mothurOutEndLine(); return 0;
		}
		
		//map sequences to bin number in the "little" otu
		map<string, int> littleBins;
        vector<string> binLabels0 = lists[0].getLabels();
		for (int i = 0; i < lists[0].getNumBins(); i++) {
		
			if (m->control_pressed) {  return 0; }
			string bin = lists[0].get(i);
            vector<string> names; m->splitAtComma(bin, names);
			for (int j = 0; j < names.size(); j++) { littleBins[names[j]] = i; }
        }
		
		ofstream out;
        map<string, string> variables; 
        variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(listFile));
        variables["[distance1]"] = lists[0].getLabel();
        variables["[tag]"] = "-"; 
        variables["[distance2]"] = lists[1].getLabel();
		string outputFileName = getOutputFileName("otuheirarchy",variables);
		m->openOutputFile(outputFileName, out);
		
		//go through each bin in "big" otu and output the bins in "little" otu which created it
        vector<string> binLabels1 = lists[1].getLabels();
		for (int i = 0; i < lists[1].getNumBins(); i++) {
		
			if (m->control_pressed) { outputTypes.clear(); out.close(); m->mothurRemove(outputFileName); return 0; }
			
			string binnames = lists[1].get(i);
            vector<string> names; m->splitAtComma(binnames, names);
			
			//output column 1
			if (output == "name")	{   out << binnames << '\t';	}
			else					{	out << binLabels1[i] << '\t';		}
			
			map<int, int> bins; //bin numbers in little that are in this bin in big
			map<int, int>::iterator it;
			
			//parse bin
			for (int j = 0; j < names.size(); j++) { bins[littleBins[names[j]]] = littleBins[names[j]];   }
			
			string col2 = "";
			for (it = bins.begin(); it != bins.end(); it++) {
				if (output == "name")	{   col2 += lists[0].get(it->first) + "\t";	}
				else					{	col2 += binLabels0[it->first] + "\t";		}
			}
			
			//output column 2
			out << col2 << endl;
		}
		
		out.close();
		
		if (m->control_pressed) { outputTypes.clear(); m->mothurRemove(outputFileName); return 0; }
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		m->mothurOut(outputFileName); m->mothurOutEndLine();	outputNames.push_back(outputFileName); outputTypes["otuheirarchy"].push_back(outputFileName); 
		m->mothurOutEndLine();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "OtuHierarchyCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************
//returns a vector of listVectors where "little" vector is first
vector<ListVector> OtuHierarchyCommand::getListVectors() {
	try {
		
		int pos; //to use in smart distancing, position of last read in file
		int lastPos;
		vector<ListVector> lists;
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;

		//open file
		ifstream in;
		m->openInputFile(listFile, in);
		
		//get first list vector in file
		ListVector* list = NULL;
		string lastLabel = "";
		if (!in.eof())	{
			pos = in.tellg();
			lastPos = pos;
			list = new ListVector(in);  
			m->gobble(in);
			lastLabel = list->getLabel();
		}
		
		while ((list != NULL) && (userLabels.size() != 0)) {
		
			if (m->control_pressed) {  in.close(); delete list;  return lists; }
			
			//is this a listvector that we want?
			if(labels.count(list->getLabel()) == 1){
				
				//make copy of listvector
				ListVector temp(*list);
				lists.push_back(temp);
			
				processedLabels.insert(list->getLabel());
				userLabels.erase(list->getLabel());
			}
		
			//you have a label the user want that is smaller than this label and the last label has not already been processed 
			if ((m->anyLabelsToProcess(list->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
				string saveLabel = list->getLabel();
				int savePos = in.tellg();
				
				//get smart distance line
				delete list;
				in.seekg(lastPos);
				if (!in.eof())	{	
					list = new ListVector(in);  
				}else { list = NULL; }
				
				//make copy of listvector
				ListVector temp(*list);
				lists.push_back(temp);
					
				processedLabels.insert(list->getLabel());
				userLabels.erase(list->getLabel());					
										
				//restore real lastlabel to save below
				list->setLabel(saveLabel);
				in.seekg(savePos);
			}
			
			lastLabel = list->getLabel();
			lastPos = pos;
			
			//get next line
			delete list;
			if (!in.eof())	{	
				pos = in.tellg();
				list = new ListVector(in);  
				m->gobble(in);
			}else { list = NULL; }
		}
		
		if (m->control_pressed) { in.close();  return lists; }				
		
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
		
		if (m->control_pressed) {  in.close(); return lists; }
		
		//run last label if you need to
		if (needToRun == true)  {
			if (list != NULL) {	delete list;	}
			
			in.seekg(lastPos);
			if (!in.eof())	{	
				list = new ListVector(in); 
				
				//make copy of listvector
				ListVector temp(*list);
				lists.push_back(temp);
				
				delete list;
			}
		}
		
		in.close();
		return lists;
	}
	catch(exception& e) {
		m->errorOut(e, "OtuHierarchyCommand", "getListVectors");
		exit(1);
	}
}

//**********************************************************************************************************************





