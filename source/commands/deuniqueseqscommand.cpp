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
        
        abort = false; calledHelp = false;
      
        vector<string> tempOutNames;
        outputTypes["fasta"] = tempOutNames;
        outputTypes["group"] = tempOutNames;
		
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
		;
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
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "DeUniqueSeqsCommand", "getOutputPattern");
        exit(1);
    }
}
/**************************************************************************************/
DeUniqueSeqsCommand::DeUniqueSeqsCommand(string option) : Command()  {	
	try {
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }
		
		else {
			OptionParser parser(option, setParameters());
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			fastaFile = validParameter.validFile(parameters, "fasta");
			if (fastaFile == "not open") { abort = true; }
			else if (fastaFile == "not found") { 				
				fastaFile = current->getFastaFile(); 
				if (fastaFile != "") { m->mothurOut("Using " + fastaFile + " as input file for the fasta parameter.\n");  }
				else { 	m->mothurOut("You have no current fastafile and the fasta parameter is required.\n");  abort = true; }
			}else { current->setFastaFile(fastaFile); }
			
			nameFile = validParameter.validFile(parameters, "name");
			if (nameFile == "not open") { abort = true; }
			else if (nameFile == "not found"){	nameFile = ""; }
			else { current->setNameFile(nameFile); }
            
            countfile = validParameter.validFile(parameters, "count");
			if (countfile == "not open") { abort = true;  }
			else if (countfile == "not found") { countfile = ""; }
			else { current->setCountFile(countfile); }
			
            if ((countfile != "") && (nameFile != "")) { m->mothurOut("When executing a deunique.seqs command you must enter ONLY ONE of the following: count or name.\n");  abort = true; }
			
            
			if ((countfile == "") && (nameFile == "")) { //look for currents
                nameFile = current->getNameFile();
				if (nameFile != "") { m->mothurOut("Using " + nameFile + " as input file for the name parameter.\n");  }
				else {
                    countfile = current->getCountFile();
                    if (countfile != "") { m->mothurOut("Using " + countfile + " as input file for the count parameter.\n");  }
                    else {  m->mothurOut("[ERROR]: You have no current name or count files one is required.\n");  abort = true; }
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
		
		if (abort) { if (calledHelp) { return 0; }  return 2;	}

		//prepare filenames and open files
		ofstream out;
        string thisOutputDir = outputdir;
		if (outputdir == "") {  thisOutputDir += util.hasPath(fastaFile);  }
        string outFastaFile = thisOutputDir + util.getRootName(util.getSimpleName(fastaFile));
       
        map<string, string> variables;
        variables["[filename]"] = outFastaFile;
        outFastaFile = getOutputFileName("fasta", variables);
		util.openOutputFile(outFastaFile, out);
		
		map<string, string> nameMap;
        CountTable ct;
        ofstream outGroup;
        string outGroupFile;
        vector<string> groups;
        if (nameFile != "") { util.readNames(nameFile, nameMap); }
        else {
            ct.readTable(countfile, true, false);
            
            if (ct.hasGroupInfo()) {
                thisOutputDir = outputdir;
                if (outputdir == "") {  thisOutputDir += util.hasPath(countfile);  }
                outGroupFile = thisOutputDir + util.getRootName(util.getSimpleName(countfile));
                variables["[filename]"] = outGroupFile;
                outGroupFile = getOutputFileName("group", variables);
                util.openOutputFile(outGroupFile, outGroup);
                groups = ct.getNamesOfGroups();
            }
        }
        
		if (m->getControl_pressed()) {  out.close(); outputTypes.clear(); util.mothurRemove(outFastaFile); if (countfile != "") { if (ct.hasGroupInfo()) { outGroup.close(); util.mothurRemove(outGroupFile); } } return 0; }
		
		ifstream in;
		util.openInputFile(fastaFile, in);
		
		while (!in.eof()) {
		
			if (m->getControl_pressed()) { in.close(); out.close(); outputTypes.clear(); util.mothurRemove(outFastaFile); if (countfile != "") { if (ct.hasGroupInfo()) { outGroup.close(); util.mothurRemove(outGroupFile); } } return 0; }
			
			Sequence seq(in); util.gobble(in);
			
			if (seq.getName() != "") {
				
                if (nameFile != "") {
                    //look for sequence name in nameMap
                    map<string, string>::iterator it = nameMap.find(seq.getName());
                    
                    if (it == nameMap.end()) {	m->mothurOut("[ERROR]: Your namefile does not contain " + seq.getName() + ", aborting.\n");  m->setControl_pressed(true); }
                    else {
                        vector<string> names;
                        util.splitAtComma(it->second, names);
                        
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
		
						
		if (m->getControl_pressed()) { outputTypes.clear(); util.mothurRemove(outFastaFile); if (countfile != "") { if (ct.hasGroupInfo()) {  util.mothurRemove(outGroupFile); } }return 0; }
		
        outputNames.push_back(outFastaFile);  outputTypes["fasta"].push_back(outFastaFile);
        if (countfile != "") { if (ct.hasGroupInfo()) { outputNames.push_back(outGroupFile);  outputTypes["group"].push_back(outGroupFile); } }
        
		m->mothurOut("\nOutput File Names: \n"); 
		for(int i = 0; i < outputNames.size(); i++) {  m->mothurOut(outputNames[i]); m->mothurOutEndLine();	 }
        m->mothurOutEndLine();
		
		
		//set fasta file as new current fastafile
		string currentName = "";
		itTypes = outputTypes.find("fasta");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setFastaFile(currentName); }
		}
        
        itTypes = outputTypes.find("group");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setGroupFile(currentName); }
		}

		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "DeUniqueSeqsCommand", "execute");
		exit(1);
	}
}
/**************************************************************************************/
