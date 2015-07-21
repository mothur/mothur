/*
 *  deuniqueseqscommand.cpp
 *  Mothur
 *
 *  Created by westcott on 10/19/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "deuniqueseqscommand.h"
#include "sequence.hpp"
#include "counttable.h"

//**********************************************************************************************************************
vector<string> DeUniqueSeqsCommand::setParameters(){	
	try {
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none","fasta",false,true,true); parameters.push_back(pfasta);
		CommandParameter pname("name", "InputTypes", "", "", "namecount", "namecount", "none","name",false,false,true); parameters.push_back(pname);
        CommandParameter pcount("count", "InputTypes", "", "", "namecount", "namecount", "none","group",false,false,true); parameters.push_back(pcount);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "DeUniqueSeqsCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string DeUniqueSeqsCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The deunique.seqs command reads a fastafile and namefile or countfile, and creates a fastafile containing all the sequences. It you provide a count file with group information a group file is also created.\n";
		helpString += "The deunique.seqs command parameters are fasta, name and count. Fasta is required and you must provide either a name or count file.\n";
		helpString += "The deunique.seqs command should be in the following format: \n";
		helpString += "deunique.seqs(fasta=yourFastaFile, name=yourNameFile) \n";	
		helpString += "Example deunique.seqs(fasta=abrecovery.unique.fasta, name=abrecovery.names).\n";
		helpString += "Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFasta).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "DeUniqueSeqsCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string DeUniqueSeqsCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "fasta") {  pattern = "[filename],redundant.fasta"; }
        else if (type == "group") {  pattern = "[filename],redundant.groups"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "DeUniqueSeqsCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
DeUniqueSeqsCommand::DeUniqueSeqsCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["fasta"] = tempOutNames;
        outputTypes["group"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "DeUniqueSeqsCommand", "DeconvoluteCommand");
		exit(1);
	}
}
/**************************************************************************************/
DeUniqueSeqsCommand::DeUniqueSeqsCommand(string option)  {	
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
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["fasta"] = tempOutNames;
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
                
                it = parameters.find("count");
				//user has given a template file
				if(it != parameters.end()){
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["count"] = inputDir + it->second;		}
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
			}
			
			nameFile = validParameter.validFile(parameters, "name", true);
			if (nameFile == "not open") { abort = true; }
			else if (nameFile == "not found"){	nameFile = ""; }
			else { m->setNameFile(nameFile); }
            
            countfile = validParameter.validFile(parameters, "count", true);
			if (countfile == "not open") { abort = true;  }
			else if (countfile == "not found") { countfile = ""; }
			else { m->setCountTableFile(countfile); }
			
            if ((countfile != "") && (nameFile != "")) { m->mothurOut("When executing a deunique.seqs command you must enter ONLY ONE of the following: count or name."); m->mothurOutEndLine(); abort = true; }
			
            
			if ((countfile == "") && (nameFile == "")) { //look for currents
                nameFile = m->getNameFile();
				if (nameFile != "") { m->mothurOut("Using " + nameFile + " as input file for the name parameter."); m->mothurOutEndLine(); }
				else {
                    countfile = m->getCountTableFile();
                    if (countfile != "") { m->mothurOut("Using " + countfile + " as input file for the count parameter."); m->mothurOutEndLine(); }
                    else {  m->mothurOut("[ERROR]: You have no current name or count files one is required."); m->mothurOutEndLine(); abort = true; }
                }
            }

		}

	}
	catch(exception& e) {
		m->errorOut(e, "DeUniqueSeqsCommand", "DeUniqueSeqsCommand");
		exit(1);
	}
}
/**************************************************************************************/
int DeUniqueSeqsCommand::execute() {	
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}

		//prepare filenames and open files
		ofstream out;
        string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += m->hasPath(fastaFile);  }
        string outFastaFile = thisOutputDir + m->getRootName(m->getSimpleName(fastaFile));
       
        map<string, string> variables;
        variables["[filename]"] = outFastaFile;
        outFastaFile = getOutputFileName("fasta", variables);
		m->openOutputFile(outFastaFile, out);
		
		map<string, string> nameMap;
        CountTable ct;
        ofstream outGroup;
        string outGroupFile;
        vector<string> groups;
        if (nameFile != "") { m->readNames(nameFile, nameMap); }
        else {
            ct.readTable(countfile, true, false);
            
            if (ct.hasGroupInfo()) {
                thisOutputDir = outputDir;
                if (outputDir == "") {  thisOutputDir += m->hasPath(countfile);  }
                outGroupFile = thisOutputDir + m->getRootName(m->getSimpleName(countfile));
                variables["[filename]"] = outGroupFile;
                outGroupFile = getOutputFileName("group", variables);
                m->openOutputFile(outGroupFile, outGroup);
                groups = ct.getNamesOfGroups();
            }
        }
        
		if (m->control_pressed) {  out.close(); outputTypes.clear(); m->mothurRemove(outFastaFile); if (countfile != "") { if (ct.hasGroupInfo()) { outGroup.close(); m->mothurRemove(outGroupFile); } } return 0; }
		
		ifstream in;
		m->openInputFile(fastaFile, in);
		
		while (!in.eof()) {
		
			if (m->control_pressed) { in.close(); out.close(); outputTypes.clear(); m->mothurRemove(outFastaFile); if (countfile != "") { if (ct.hasGroupInfo()) { outGroup.close(); m->mothurRemove(outGroupFile); } } return 0; }
			
			Sequence seq(in); m->gobble(in);
			
			if (seq.getName() != "") {
				
                if (nameFile != "") {
                    //look for sequence name in nameMap
                    map<string, string>::iterator it = nameMap.find(seq.getName());
                    
                    if (it == nameMap.end()) {	m->mothurOut("[ERROR]: Your namefile does not contain " + seq.getName() + ", aborting."); m->mothurOutEndLine(); m->control_pressed = true; }
                    else {
                        vector<string> names;
                        m->splitAtComma(it->second, names);
                        
                        //output sequences
                        for (int i = 0; i < names.size(); i++) {
                            out << ">" << names[i] << endl;
                            out << seq.getAligned() << endl;
                        }
                        
                        //remove seq from name map so we can check for seqs in namefile not in fastafile later
                        nameMap.erase(it);
                    }
                }else {
                    if (ct.hasGroupInfo()) {
                        vector<int> groupCounts = ct.getGroupCounts(seq.getName());
                        int count = 1;
                        for (int i = 0; i < groups.size(); i++) {
                            for (int j = 0; j < groupCounts[i]; j++) {
                                outGroup << seq.getName()+"_"+toString(count) << '\t' << groups[i] << endl; count++;
                            }
                        }
                        
                    }
                    
                    int numReps = ct.getNumSeqs(seq.getName()); //will report error and set m->control_pressed if not found
                    for (int i = 0; i < numReps; i++) {
                        out << ">" << seq.getName()+"_"+toString(i+1) << endl;
                        out << seq.getAligned() << endl;
                    }
                }
			}
		}
		in.close();
		out.close();
        if (countfile != "") { if (ct.hasGroupInfo()) { outGroup.close(); } }
		
						
		if (m->control_pressed) { outputTypes.clear(); m->mothurRemove(outFastaFile); if (countfile != "") { if (ct.hasGroupInfo()) {  m->mothurRemove(outGroupFile); } }return 0; }
		
        outputNames.push_back(outFastaFile);  outputTypes["fasta"].push_back(outFastaFile);
        if (countfile != "") { if (ct.hasGroupInfo()) { outputNames.push_back(outGroupFile);  outputTypes["group"].push_back(outGroupFile); } }
        
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for(int i = 0; i < outputNames.size(); i++) {  m->mothurOut(outputNames[i]); m->mothurOutEndLine();	 }
        m->mothurOutEndLine();
		
		
		//set fasta file as new current fastafile
		string current = "";
		itTypes = outputTypes.find("fasta");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setFastaFile(current); }
		}
        
        itTypes = outputTypes.find("group");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setGroupFile(current); }
		}

		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "DeUniqueSeqsCommand", "execute");
		exit(1);
	}
}
/**************************************************************************************/
