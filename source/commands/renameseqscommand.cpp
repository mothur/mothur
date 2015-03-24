//
//  renameseqscommand.cpp
//  Mothur
//
//  Created by SarahsWork on 5/28/13.
//  Copyright (c) 2013 Schloss Lab. All rights reserved.
//

#include "renameseqscommand.h"
#include "sequence.hpp"
#include "groupmap.h"

//**********************************************************************************************************************
vector<string> RenameSeqsCommand::setParameters(){
	try {
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none","fasta",false,true,true); parameters.push_back(pfasta);
        CommandParameter pname("name", "InputTypes", "", "", "NameCount", "none", "none","name",false,false,true); parameters.push_back(pname);
		CommandParameter pgroup("group", "InputTypes", "", "", "CountGroup", "none", "none","group",false,true,true); parameters.push_back(pgroup);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "RenameSeqsCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string RenameSeqsCommand::getHelpString(){
	try {
		string helpString = "";
		helpString += "The rename.seqs command reads a fastafile and groupfile with an optional namefile, and creates files with the sequence names concatenated with the group. For example if a line in the group file is 'seq1   group1', the new sequence name will be seq1_group1.\n";
		helpString += "The rename.seqs command parameters are fasta, name and group. Fasta and group are required, unless a current file is available for both.\n";
		helpString += "The rename.seqs command should be in the following format: \n";
		helpString += "rename.seqs(fasta=yourFastaFile, group=yourGroupFile) \n";
		helpString += "Example rename.seqs(fasta=abrecovery.unique.fasta, group=abrecovery.group).\n";
		helpString += "Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFasta).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "RenameSeqsCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string RenameSeqsCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "fasta")        {  pattern = "[filename],renamed,[extension]"; }
        else if (type == "name")    {  pattern = "[filename],renamed,[extension]"; }
        else if (type == "group")   {  pattern = "[filename],renamed,[extension]"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "RenameSeqsCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
RenameSeqsCommand::RenameSeqsCommand(){
	try {
		abort = true; calledHelp = true;
		setParameters();
		vector<string> tempOutNames;
		outputTypes["fasta"] = tempOutNames;
        outputTypes["name"] = tempOutNames;
        outputTypes["group"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "RenameSeqsCommand", "RenameSeqsCommand");
		exit(1);
	}
}
/**************************************************************************************/
RenameSeqsCommand::RenameSeqsCommand(string option)  {
	try {
		abort = false; calledHelp = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			map<string, string>::iterator it;
            
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) {
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			vector<string> tempOutNames;
            outputTypes["fasta"] = tempOutNames;
            outputTypes["name"] = tempOutNames;
            outputTypes["group"] = tempOutNames;
            
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

			}
            
			
			//check for required parameters
			fastaFile = validParameter.validFile(parameters, "fasta", true);
			if (fastaFile == "not open") { abort = true; }
			else if (fastaFile == "not found") {
				fastaFile = m->getFastaFile();
				if (fastaFile != "") { m->mothurOut("Using " + fastaFile + " as input file for the fasta parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current fastafile and the fasta parameter is required."); m->mothurOutEndLine(); abort = true; }
			}else { m->setFastaFile(fastaFile); }
			
			//if the user changes the output directory command factory will send this info to us in the output parameter
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){
				outputDir = "";
				outputDir += m->hasPath(fastaFile); //if user entered a file with a path then preserve it
			}
			
            groupfile = validParameter.validFile(parameters, "group", true);
            if (groupfile == "not open") { abort = true; }
			else if (groupfile == "not found") {
				groupfile = m->getGroupFile();
				if (groupfile != "") { m->mothurOut("Using " + groupfile + " as input file for the group parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current groupfile and the group parameter is required."); m->mothurOutEndLine(); abort = true; }
			}else { m->setGroupFile(groupfile); }
			
			//if the user changes the output directory command factory will send this info to us in the output parameter
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){ outputDir = ""; }

			nameFile = validParameter.validFile(parameters, "name", true);
			if (nameFile == "not open") { abort = true; }
			else if (nameFile == "not found"){ nameFile =""; }
			else { m->setNameFile(nameFile); }
            
            if (nameFile == "") {
                vector<string> files; files.push_back(fastaFile);
                parser.getNameFile(files);
            }

		}
        
	}
	catch(exception& e) {
		m->errorOut(e, "RenameSeqsCommand", "RenameSeqsCommand");
		exit(1);
	}
}
/**************************************************************************************/
int RenameSeqsCommand::execute() {
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
        
        GroupMap groupMap(groupfile);
        groupMap.readMap();
                
		//prepare filenames and open files
        string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += m->hasPath(fastaFile);  }
		string outFastaFile = thisOutputDir + m->getRootName(m->getSimpleName(fastaFile));
        map<string, string> variables;
        variables["[filename]"] = outFastaFile;
        variables["[extension]"] = m->getExtension(fastaFile);
        outFastaFile = getOutputFileName("fasta", variables);
        outputNames.push_back(outFastaFile); outputTypes["fasta"].push_back(outFastaFile);
        
        ofstream outFasta;
		m->openOutputFile(outFastaFile, outFasta);
        
        ifstream in;
        m->openInputFile(fastaFile, in);
        
        while (!in.eof()) {
            if (m->control_pressed) { break; }
            
            Sequence seq(in); m->gobble(in);
            string group = groupMap.getGroup(seq.getName());
            if (group == "not found") {  m->mothurOut("[ERROR]: " + seq.getName() + " is not in your group file, please correct.\n"); m->control_pressed = true; }
            else {
                string newName = seq.getName() + "_" + group;
                seq.setName(newName);
                seq.printSequence(outFasta);
            }
            
        }
        in.close();
        
        if (m->control_pressed) {  for (int i = 0; i < outputNames.size(); i++) { m->mothurRemove(outputNames[i]);  } return 0; }
        
        bool notDone = true;
        if (nameFile != "") {
            thisOutputDir = outputDir;
            if (outputDir == "") {  thisOutputDir += m->hasPath(nameFile);  }
            string outNameFile = thisOutputDir + m->getRootName(m->getSimpleName(nameFile));
            variables["[filename]"] = outNameFile;
            variables["[extension]"] = m->getExtension(nameFile);
            outNameFile = getOutputFileName("group", variables);
            outputNames.push_back(outNameFile); outputTypes["name"].push_back(outNameFile);
            
            ofstream outName;
            m->openOutputFile(outNameFile, outName);
            
            map<string, vector<string> > nameMap;
            m->readNames(nameFile, nameMap);
            
            //process name file changing names
            for (map<string, vector<string> >::iterator it = nameMap.begin(); it != nameMap.end(); it++) {
                for (int i = 0; i < (it->second).size()-1; i++) {
                    if (m->control_pressed) { break; }
                    string group = groupMap.getGroup((it->second)[i]);
                    if (group == "not found") {  m->mothurOut("[ERROR]: " + (it->second)[i] + " is not in your group file, please correct.\n"); m->control_pressed = true;  }
                    else {
                        string newName = (it->second)[i] + "_" + group;
                        groupMap.renameSeq((it->second)[i], newName); //change in group file
                        (it->second)[i] = newName; //change in namefile
                    }
                    if (i == 0) {  outName << (it->second)[i] << '\t' << (it->second)[i] << ','; }
                    else { outName << (it->second)[i] << ','; }
                }
                
                //print last one
                if ((it->second).size() == 1) {
                    string group = groupMap.getGroup((it->second)[0]);
                    if (group == "not found") {  m->mothurOut("[ERROR]: " + (it->second)[0] + " is not in your group file, please correct.\n"); m->control_pressed = true;  }
                    else {
                        string newName = (it->second)[0] + "_" + group;
                        groupMap.renameSeq((it->second)[0], newName); //change in group file
                        (it->second)[0] = newName; //change in namefile

                        outName << (it->second)[0] << '\t' << (it->second)[0] << endl;
                    }
                }
                else {
                    string group = groupMap.getGroup((it->second)[(it->second).size()-1]);
                    if (group == "not found") {  m->mothurOut("[ERROR]: " + (it->second)[(it->second).size()-1] + " is not in your group file, please correct.\n"); m->control_pressed = true;  }
                    else {
                        string newName = (it->second)[(it->second).size()-1] + "_" + group;
                        groupMap.renameSeq((it->second)[(it->second).size()-1], newName); //change in group file
                        (it->second)[(it->second).size()-1] = newName; //change in namefile

                        outName << (it->second)[(it->second).size()-1] << endl;
                    }
                }
            }
            notDone = false;
            outName.close();
        }
        
        if (m->control_pressed) {  for (int i = 0; i < outputNames.size(); i++) { m->mothurRemove(outputNames[i]);  } return 0; }
        
        if (notDone) {
            vector<string> seqs = groupMap.getNamesSeqs();
            for (int i = 0; i < seqs.size(); i++) {
                if (m->control_pressed) { break; }
                string group = groupMap.getGroup(seqs[i]);
                string newName = seqs[i] + "_" + group;
                groupMap.renameSeq(seqs[i], newName);
            }
        }
        if (m->control_pressed) {  for (int i = 0; i < outputNames.size(); i++) { m->mothurRemove(outputNames[i]);  } return 0; }
        
        thisOutputDir = outputDir;
        if (outputDir == "") {  thisOutputDir += m->hasPath(groupfile);  }
		string outGroupFile = thisOutputDir + m->getRootName(m->getSimpleName(groupfile));
        variables["[filename]"] = outGroupFile;
        variables["[extension]"] = m->getExtension(groupfile);
        outGroupFile = getOutputFileName("group", variables);
        outputNames.push_back(outGroupFile); outputTypes["group"].push_back(outGroupFile);
        
        ofstream outGroup;
		m->openOutputFile(outGroupFile, outGroup);
        groupMap.print(outGroup);
        outGroup.close();
        
        if (m->control_pressed) {  for (int i = 0; i < outputNames.size(); i++) { m->mothurRemove(outputNames[i]);  } return 0; }

        m->mothurOutEndLine();
        m->mothurOut("Output File Names: "); m->mothurOutEndLine();
        for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
        m->mothurOutEndLine();
        
        //set fasta file as new current fastafile
        string current = "";
        itTypes = outputTypes.find("fasta");
        if (itTypes != outputTypes.end()) {
            if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setFastaFile(current); }
        }
        
        itTypes = outputTypes.find("name");
        if (itTypes != outputTypes.end()) {
            if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setNameFile(current); }
        }
        
        itTypes = outputTypes.find("group");
        if (itTypes != outputTypes.end()) {
            if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setGroupFile(current); }
        }
				
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "RenameSeqsCommand", "execute");
		exit(1);
	}
}
/**************************************************************************************/

