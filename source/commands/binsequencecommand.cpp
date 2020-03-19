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
        else if (type == "count") {  pattern = "[filename],count_table";  }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
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
        outputTypes["count"] = tempOutNames;
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
		allLines = true;
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
				if (!validParameter.isValidParameter(it->first, myArray, it->second)) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["fasta"] = tempOutNames;
            outputTypes["count"] = tempOutNames;
			
			//check for required parameters
			fastafile = validParameter.validFile(parameters, "fasta");
			if (fastafile == "not found") { 				//if there is a current phylip file, use it
				fastafile = current->getFastaFile(); 
				if (fastafile != "") { m->mothurOut("Using " + fastafile + " as input file for the fasta parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current fasta file and the fasta parameter is required."); m->mothurOutEndLine(); abort = true; }
			}
			else if (fastafile == "not open") { abort = true; }	
			else { current->setFastaFile(fastafile); }
			
			listfile = validParameter.validFile(parameters, "list");
			if (listfile == "not found") { 			
				listfile = current->getListFile(); 
				if (listfile != "") { m->mothurOut("Using " + listfile + " as input file for the list parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current list file and the list parameter is required."); m->mothurOutEndLine(); abort = true; }
			}
			else if (listfile == "not open") { listfile = ""; abort = true; }	
			else { current->setListFile(listfile); }
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){
				outputDir = "";	
				outputDir += util.hasPath(listfile); //if user entered a file with a path then preserve it	
			}
			
		
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			
			label = validParameter.valid(parameters, "label");			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  util.splitAtDash(label, labels);  allLines = false;  }
				else { allLines = true;  }
			}
			
			string namesfile = validParameter.validFile(parameters, "name");
			if (namesfile == "not open") { namesfile = ""; abort = true; }	
			else if (namesfile == "not found") { namesfile = ""; }
			else {  current->setNameFile(namesfile); }

			string groupfile = validParameter.validFile(parameters, "group");
			if (groupfile == "not open") { abort = true; }
			else if (groupfile == "not found") { groupfile = ""; }
			else { current->setGroupFile(groupfile); }
            
            countfile = validParameter.validFile(parameters, "count");
			if (countfile == "not open") { countfile = ""; abort = true; }
			else if (countfile == "not found") { countfile = "";  }	
			else { current->setCountFile(countfile); }
            
            if ((namesfile != "") && (countfile != "")) { m->mothurOut("[ERROR]: you may only use one of the following: name or count.\n");  abort = true; }
            if ((groupfile != "") && (countfile != "")) {  m->mothurOut("[ERROR]: you may only use one of the following: group or count.\n");  abort=true; }
            
            if (!abort) {
            if ((namesfile != "") || (groupfile != "")) { //convert to count
                
                string rootFileName = namesfile;
                if (rootFileName == "") { rootFileName = groupfile; }
                
                if (outputDir == "") { outputDir = util.hasPath(rootFileName); }
                map<string, string> variables; variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(rootFileName));
                string outputFileName = getOutputFileName("count", variables);

                CountTable ct; ct.createTable(namesfile, groupfile, nullVector); ct.printCompressedTable(outputFileName);
                outputNames.push_back(outputFileName); outputTypes["count"].push_back(outputFileName);
                
                current->setCountFile(outputFileName);
                countfile = outputFileName;
                
                //list file will contain redund names since name file is provided - remove dups
                string tempAccnos = namesfile + ".accnos.temp";
                vector<string> namesOfSeqs = ct.getNamesOfSeqs();
                util.printAccnos(tempAccnos, namesOfSeqs);
                
                string inputString = "list=" + listfile + ", accnos=" + tempAccnos;
                
                m->mothurOut("/******************************************/\n");
                m->mothurOut("\nRunning command: get.seqs(" + inputString + ")\n");
                current->setMothurCalling(true);
                
                Command* getSeqsCommand = new GetSeqsCommand(inputString);
                getSeqsCommand->execute();
                
                string templistfile = getSeqsCommand->getOutputFiles()["list"][0];
                string newName = util.getRootName(listfile) + "unique.list";
                util.renameFile(templistfile, newName);  listfile = newName;
                namesfile = ""; groupfile = ""; util.mothurRemove(tempAccnos);
                current->setListFile(listfile);
                
                delete getSeqsCommand;
                current->setMothurCalling(false);
                
                m->mothurOut("/******************************************/\n");
            }
            }
			
            if (countfile == "") {
                if (namesfile == ""){
                    vector<string> files; files.push_back(fastafile); 
                    if (!current->getMothurCalling())  {  parser.getNameFile(files);  }
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
		if (abort) { if (calledHelp) { return 0; }  return 2;	}

		FastaMap fasta; fasta.readFastaFile(fastafile);
        
		//if user gave a namesfile then use it
        if (countfile != "") {  ct.readTable(countfile, true, false);  }
		
		InputData input(listfile, "list", nullVector);
        set<string> processedLabels;
        set<string> userLabels = labels;
        string lastLabel = "";
        int error = 0;
        
        ListVector* list = util.getNextList(input, allLines, userLabels, processedLabels, lastLabel);
        
        while (list != NULL) {
            
            if (m->getControl_pressed()) { delete list; break; }
            
            error = process(list, fasta); delete list;
            if (error == 1)  { for (int i = 0; i < outputNames.size(); i++) {    util.mothurRemove(outputNames[i]);        }  return 0; }
            
            list = util.getNextList(input, allLines, userLabels, processedLabels, lastLabel);
        }

		if(m->getControl_pressed())  { for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]);		}  return 0; }
        
        //set align file as new current fastafile
		string currentFasta = "";
		itTypes = outputTypes.find("fasta");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentFasta = (itTypes->second)[0]; current->setFastaFile(currentFasta); }
		}
		
		m->mothurOut("\nOutput File Names: \n"); 
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i] +"\n"); 	} m->mothurOutEndLine();

		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "BinSeqCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
//return 1 if error, 0 otherwise
int BinSeqCommand::process(ListVector* list, FastaMap& fasta) {
	try {
        map<string, string> variables; 
        variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(listfile));
        variables["[distance]"] = list->getLabel();
        string outputFileName = getOutputFileName("fasta", variables);
        
        ofstream out; util.openOutputFile(outputFileName, out);
        outputNames.push_back(outputFileName);  outputTypes["fasta"].push_back(outputFileName);
        
        m->mothurOut(list->getLabel()); m->mothurOutEndLine();
        
        //for each bin in the list vector
        vector<string> binLabels = list->getLabels();
        for (int i = 0; i < list->size(); i++) {
            
            if (m->getControl_pressed()) {  return 1; }
            
            string binnames = list->get(i);
            vector<string> names;
            util.splitAtComma(binnames, names);
            
            for (int j = 0; j < names.size(); j++) {
                string name = names[j];
                
                //do work for that name
                string sequence = fasta.getSequence(name);
                
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
                        
                    }else { m->mothurOut(name + " is missing from your fasta. Does your list file contain all sequence names or just the uniques?\n"); return 1; }
                }else {
                    if (sequence != "not found") {
                       
                        name = name + "\t" + binLabels[i];
                        out << ">" << name << endl;
                        out << sequence << endl;
                        
                    }else { m->mothurOut(name + " is missing from your fasta or count file. Please correct. \n");  return 1; }
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


