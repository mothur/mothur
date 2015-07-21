/*
 *  binsequencecommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 4/3/09.
 *  Copyright 2009 Schloss Lab UMASS Amhers. All rights reserved.
 *
 */

#include "binsequencecommand.h"


//**********************************************************************************************************************
vector<string> BinSeqCommand::setParameters(){	
	try {
		CommandParameter plist("list", "InputTypes", "", "", "none", "none", "none","",false,true,true); parameters.push_back(plist);
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none","fasta",false,true,true); parameters.push_back(pfasta);
        CommandParameter pname("name", "InputTypes", "", "", "NameCount", "none", "none","",false,false,true); parameters.push_back(pname);
        CommandParameter pcount("count", "InputTypes", "", "", "NameCount-CountGroup", "none", "none","",false,false,true); parameters.push_back(pcount);
		CommandParameter pgroup("group", "InputTypes", "", "", "CountGroup", "none", "none","",false,false,true); parameters.push_back(pgroup);
		CommandParameter plabel("label", "String", "", "", "", "", "","",false,false); parameters.push_back(plabel);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
	
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "BinSeqCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string BinSeqCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The bin.seqs command parameters are list, fasta, name, count, label and group.  The fasta and list are required, unless you have a valid current list and fasta file.\n";
		helpString += "The label parameter allows you to select what distance levels you would like a output files created for, and are separated by dashes.\n";
		helpString += "The bin.seqs command should be in the following format: bin.seqs(fasta=yourFastaFile, name=yourNamesFile, group=yourGroupFile, label=yourLabels).\n";
		helpString += "Example bin.seqs(fasta=amazon.fasta, group=amazon.groups, name=amazon.names).\n";
		helpString += "The default value for label is all lines in your inputfile.\n";
		helpString += "The bin.seqs command outputs a .fasta file for each distance you specify appending the OTU number to each name.\n";
		helpString += "If you provide a groupfile, then it also appends the sequences group to the name.\n";
		helpString += "Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFastaFile).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "BinSeqCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string BinSeqCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "fasta") {  pattern = "[filename],[distance],fasta"; } //makes file like: amazon.0.03.fasta
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "BinSeqCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
BinSeqCommand::BinSeqCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["fasta"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "BinSeqCommand", "BinSeqCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
BinSeqCommand::BinSeqCommand(string option) {
	try {
		abort = false; calledHelp = false;   
		allLines = 1;
		labels.clear();
		
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
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["fasta"] = tempOutNames;
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("fasta");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["fasta"] = inputDir + it->second;		}
				}
				
				it = parameters.find("list");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["list"] = inputDir + it->second;		}
				}
				
				it = parameters.find("name");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["name"] = inputDir + it->second;		}
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
			}

			
			//check for required parameters
			fastafile = validParameter.validFile(parameters, "fasta", true);
			if (fastafile == "not found") { 				//if there is a current phylip file, use it
				fastafile = m->getFastaFile(); 
				if (fastafile != "") { m->mothurOut("Using " + fastafile + " as input file for the fasta parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current fasta file and the fasta parameter is required."); m->mothurOutEndLine(); abort = true; }
			}
			else if (fastafile == "not open") { abort = true; }	
			else { m->setFastaFile(fastafile); }
			
			listfile = validParameter.validFile(parameters, "list", true);
			if (listfile == "not found") { 			
				listfile = m->getListFile(); 
				if (listfile != "") { m->mothurOut("Using " + listfile + " as input file for the list parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current list file and the list parameter is required."); m->mothurOutEndLine(); abort = true; }
			}
			else if (listfile == "not open") { listfile = ""; abort = true; }	
			else { m->setListFile(listfile); }
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += m->hasPath(listfile); //if user entered a file with a path then preserve it	
			}
			
		
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			
			label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  m->splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
			
			namesfile = validParameter.validFile(parameters, "name", true);
			if (namesfile == "not open") { namesfile = ""; abort = true; }	
			else if (namesfile == "not found") { namesfile = ""; }
			else {  m->setNameFile(namesfile); }

			groupfile = validParameter.validFile(parameters, "group", true);
			if (groupfile == "not open") { abort = true; }
			else if (groupfile == "not found") { groupfile = ""; }
			else { m->setGroupFile(groupfile); }
            
            countfile = validParameter.validFile(parameters, "count", true);
			if (countfile == "not open") { countfile = ""; abort = true; }
			else if (countfile == "not found") { countfile = "";  }	
			else { m->setCountTableFile(countfile); }
            
            if ((namesfile != "") && (countfile != "")) {
                m->mothurOut("[ERROR]: you may only use one of the following: name or count."); m->mothurOutEndLine(); abort = true;
            }
			
            if ((groupfile != "") && (countfile != "")) {
                m->mothurOut("[ERROR]: you may only use one of the following: group or count."); m->mothurOutEndLine(); abort=true;
            }
			
            if (countfile == "") {
                if (namesfile == ""){
                    vector<string> files; files.push_back(fastafile); 
                    parser.getNameFile(files);
                }
            }
			
		}
	}
	catch(exception& e) {
		m->errorOut(e, "BinSeqCommand", "BinSeqCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

BinSeqCommand::~BinSeqCommand(){}
//**********************************************************************************************************************

int BinSeqCommand::execute(){
	try {
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
	
		int error = 0;
		
		fasta = new FastaMap();
		if (groupfile != "") {
			groupMap = new GroupMap(groupfile);
			groupMap->readMap();
		}
		
		//read fastafile
		fasta->readFastaFile(fastafile);
		
		//if user gave a namesfile then use it
		if (namesfile != "") {  readNamesFile();  }
        if (countfile != "") {  ct.readTable(countfile, true, false);  }
		
		input = new InputData(listfile, "list");
		list = input->getListVector();
		string lastLabel = list->getLabel();
		
		if (m->control_pressed) {  delete input;  delete fasta; if (groupfile != "") {  delete groupMap;   } return 0; }
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;

				
		while((list != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			
			if(m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);		} delete input;  delete fasta; if (groupfile != "") {  delete groupMap;   } return 0; }	
			
			if(allLines == 1 || labels.count(list->getLabel()) == 1){
				
				error = process(list);	
				if (error == 1) { for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);		} delete input;  delete fasta; if (groupfile != "") {  delete groupMap;   } return 0; }	
							
				processedLabels.insert(list->getLabel());
				userLabels.erase(list->getLabel());
			}
			
			if ((m->anyLabelsToProcess(list->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
				string saveLabel = list->getLabel();
				
				delete list;
				list = input->getListVector(lastLabel);
				
				error = process(list);	
				if (error == 1) { for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);		} delete input;  delete fasta; if (groupfile != "") {  delete groupMap;   } return 0; }
													
				processedLabels.insert(list->getLabel());
				userLabels.erase(list->getLabel());
				
				//restore real lastlabel to save below
				list->setLabel(saveLabel);
			}
			
			lastLabel = list->getLabel();			
			
			delete list;
			list = input->getListVector();
		}
		
		if(m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);		} delete input;  delete fasta; if (groupfile != "") {  delete groupMap;   } return 0; }	

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
		if (needToRun == true)  {
			if (list != NULL) {		delete list;	}
			list = input->getListVector(lastLabel);
				
			error = process(list);	
			if (error == 1) { for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);		} delete input;  delete fasta; if (groupfile != "") {  delete groupMap;   } return 0; }
			
			delete list;  
		}
		
		delete input;  
		delete fasta; 
		if (groupfile != "") {  delete groupMap;   } 
		
		if(m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);		}  return 0; }	
        
        //set align file as new current fastafile
		string currentFasta = "";
		itTypes = outputTypes.find("fasta");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentFasta = (itTypes->second)[0]; m->setFastaFile(currentFasta); }
		}
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();

		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "BinSeqCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************
void BinSeqCommand::readNamesFile() {
	try {
		vector<string> dupNames;
		m->openInputFile(namesfile, inNames);
		
		string name, names, sequence;
	
		while(inNames){
			inNames >> name;			//read from first column  A
			inNames >> names;		//read from second column  A,B,C,D
			
			dupNames.clear();
			
			//parse names into vector
			m->splitAtComma(names, dupNames);
			
			//store names in fasta map
			sequence = fasta->getSequence(name);
			for (int i = 0; i < dupNames.size(); i++) {
				fasta->push_back(dupNames[i], sequence);
			}
		
			m->gobble(inNames);
		}
		inNames.close();

	}
	catch(exception& e) {
		m->errorOut(e, "BinSeqCommand", "readNamesFile");
		exit(1);
	}
}
//**********************************************************************************************************************
//return 1 if error, 0 otherwise
int BinSeqCommand::process(ListVector* list) {
	try {
        map<string, string> variables; 
        variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(listfile));
        variables["[distance]"] = list->getLabel();
        string outputFileName = getOutputFileName("fasta", variables);
        
        m->openOutputFile(outputFileName, out);
        outputNames.push_back(outputFileName);  outputTypes["fasta"].push_back(outputFileName);
        
        m->mothurOut(list->getLabel()); m->mothurOutEndLine();
        
        //for each bin in the list vector
        vector<string> binLabels = list->getLabels();
        for (int i = 0; i < list->size(); i++) {
            
            if (m->control_pressed) {  return 1; }
            
            string binnames = list->get(i);
            vector<string> names;
            m->splitAtComma(binnames, names);
            for (int j = 0; j < names.size(); j++) {
                string name = names[j];
                
                //do work for that name
                string sequence = fasta->getSequence(name);
                
                if (countfile != "") {
                    if (sequence != "not found") {
                        if (ct.hasGroupInfo()) {
                            vector<string> groups = ct.getGroups(name);
                            string groupInfo = "";
                            for (int k = 0; k < groups.size()-1; k++) {
                                groupInfo += groups[k] + "-";
                            }
                            if (groups.size() != 0) { groupInfo += groups[groups.size()-1]; }
                            else { groupInfo = "not found";  }
                            name = name + "\t" + groupInfo + "\t" + binLabels[i] + "\tNumRep=" + toString(ct.getNumSeqs(name));
                            out << ">" << name << endl;
                            out << sequence << endl;
                        }else {
                            name = name + "\t" + binLabels[i] + "\tNumRep=" + toString(ct.getNumSeqs(name));
                            out << ">" << name << endl;
                            out << sequence << endl;
                        }
                        
                    }else { m->mothurOut(name + " is missing from your fasta. Does your list file contain all sequence names or just the uniques?"); m->mothurOutEndLine(); return 1; }
                }else {
                    if (sequence != "not found") {
                        //if you don't have groups
                        if (groupfile == "") {
                            name = name + "\t" + binLabels[i];
                            out << ">" << name << endl;
                            out << sequence << endl;
                        }else {//if you do have groups
                            string group = groupMap->getGroup(name);
                            if (group == "not found") {  
                                m->mothurOut(name + " is missing from your group file. Please correct. ");  m->mothurOutEndLine();
                                return 1;
                            }else{
                                name = name + "\t" + group + "\t" + binLabels[i];
                                out << ">" << name << endl;
                                out << sequence << endl;
                            }
                        }
                    }else { m->mothurOut(name + " is missing from your fasta or name file. Please correct. "); m->mothurOutEndLine(); return 1; }
                }
            }
        }
        
        out.close();
        return 0;

	}
	catch(exception& e) {
		m->errorOut(e, "BinSeqCommand", "process");
		exit(1);
	}
}
//**********************************************************************************************************************


