/*
 *  deconvolute.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 1/21/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "deconvolutecommand.h"
#include "sequence.hpp"

//**********************************************************************************************************************
vector<string> DeconvoluteCommand::setParameters(){	
	try {
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(pfasta);
		CommandParameter pname("name", "InputTypes", "", "", "namecount", "none", "none",false,false); parameters.push_back(pname);
        CommandParameter pcount("count", "InputTypes", "", "", "namecount", "none", "none",false,false); parameters.push_back(pcount);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "DeconvoluteCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string DeconvoluteCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The unique.seqs command reads a fastafile and creates a name or count file.\n";
		helpString += "It creates a file where the first column is the groupname and the second column is a list of sequence names who have the same sequence. \n";
		helpString += "If the sequence is unique the second column will just contain its name. \n";
		helpString += "The unique.seqs command parameters are fasta and name.  fasta is required, unless there is a valid current fasta file.\n";
		helpString += "The unique.seqs command should be in the following format: \n";
		helpString += "unique.seqs(fasta=yourFastaFile) \n";	
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "DeconvoluteCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string DeconvoluteCommand::getOutputFileNameTag(string type, string inputName=""){	
	try {
        string outputFileName = "";
		map<string, vector<string> >::iterator it;
        
        //is this a type this command creates
        it = outputTypes.find(type);
        if (it == outputTypes.end()) {  m->mothurOut("[ERROR]: this command doesn't create a " + type + " output file.\n"); }
        else {
            if (type == "fasta") {  outputFileName =  "unique" + m->getExtension(inputName); }
            else if (type == "name") {  outputFileName =  "names"; }
            else if (type == "count") {  outputFileName =  "count.table"; }
            else { m->mothurOut("[ERROR]: No definition for type " + type + " output file tag.\n"); m->control_pressed = true;  }
        }
        return outputFileName;
	}
	catch(exception& e) {
		m->errorOut(e, "DeconvoluteCommand", "getOutputFileNameTag");
		exit(1);
	}
}
//**********************************************************************************************************************
DeconvoluteCommand::DeconvoluteCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["fasta"] = tempOutNames;
		outputTypes["name"] = tempOutNames;
        outputTypes["count"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "DeconvoluteCommand", "DeconvoluteCommand");
		exit(1);
	}
}
/**************************************************************************************/
DeconvoluteCommand::DeconvoluteCommand(string option)  {	
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
			outputTypes["name"] = tempOutNames;
            outputTypes["count"] = tempOutNames;
		
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
			inFastaName = validParameter.validFile(parameters, "fasta", true);
			if (inFastaName == "not open") { abort = true; }
			else if (inFastaName == "not found") { 				
				inFastaName = m->getFastaFile(); 
				if (inFastaName != "") { m->mothurOut("Using " + inFastaName + " as input file for the fasta parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current fastafile and the fasta parameter is required."); m->mothurOutEndLine(); abort = true; }
			}else { m->setFastaFile(inFastaName); }	
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += m->hasPath(inFastaName); //if user entered a file with a path then preserve it	
			}
			
			oldNameMapFName = validParameter.validFile(parameters, "name", true);
			if (oldNameMapFName == "not open") { oldNameMapFName = ""; abort = true; }
			else if (oldNameMapFName == "not found"){	oldNameMapFName = "";	}
			else { m->setNameFile(oldNameMapFName); }
            
            countfile = validParameter.validFile(parameters, "count", true);
			if (countfile == "not open") { abort = true; countfile = ""; }	
			else if (countfile == "not found") { countfile = ""; }
			else { m->setCountTableFile(countfile); }
			
            if ((countfile != "") && (oldNameMapFName != "")) { m->mothurOut("When executing a unique.seqs command you must enter ONLY ONE of the following: count or name."); m->mothurOutEndLine(); abort = true; }
			

			if (countfile == "") {
                if (oldNameMapFName == "") {
                    vector<string> files; files.push_back(inFastaName);
                    parser.getNameFile(files);
                }
            }
			
		}

	}
	catch(exception& e) {
		m->errorOut(e, "DeconvoluteCommand", "DeconvoluteCommand");
		exit(1);
	}
}
/**************************************************************************************/
int DeconvoluteCommand::execute() {	
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}

		//prepare filenames and open files
		string outNameFile = outputDir + m->getRootName(m->getSimpleName(inFastaName)) + getOutputFileNameTag("name");
        string outCountFile = outputDir + m->getRootName(m->getSimpleName(inFastaName)) + getOutputFileNameTag("count");
		string outFastaFile = outputDir + m->getRootName(m->getSimpleName(inFastaName)) + getOutputFileNameTag("fasta", inFastaName);
		
		map<string, string> nameMap;
		map<string, string>::iterator itNames;
		if (oldNameMapFName != "")  {  
            m->readNames(oldNameMapFName, nameMap); 
            if (oldNameMapFName == outNameFile){ outNameFile = outputDir + m->getRootName(m->getSimpleName(inFastaName)) + "unique." + getOutputFileNameTag("name");   }
        }
        CountTable ct;
        if (countfile != "")  {  
            ct.readTable(countfile);
            if (countfile == outCountFile){ outCountFile = outputDir + m->getRootName(m->getSimpleName(inFastaName)) + "unique." + getOutputFileNameTag("count");   }
        }
		
		if (m->control_pressed) { return 0; }
		
		ifstream in; 
		m->openInputFile(inFastaName, in);
		
		ofstream outFasta;
		m->openOutputFile(outFastaFile, outFasta);
		
		map<string, string> sequenceStrings; //sequenceString -> list of names.  "atgc...." -> seq1,seq2,seq3.
		map<string, string>::iterator itStrings;
		set<string> nameInFastaFile; //for sanity checking
		set<string>::iterator itname;
		vector<string> nameFileOrder;
		int count = 0;
		while (!in.eof()) {
			
			if (m->control_pressed) { in.close(); outFasta.close(); m->mothurRemove(outFastaFile); return 0; }
			
			Sequence seq(in);
			
			if (seq.getName() != "") {
				
				//sanity checks
				itname = nameInFastaFile.find(seq.getName());
				if (itname == nameInFastaFile.end()) { nameInFastaFile.insert(seq.getName());  }
				else { m->mothurOut("[ERROR]: You already have a sequence named " + seq.getName() + " in your fasta file, sequence names must be unique, please correct."); m->mothurOutEndLine(); }

				itStrings = sequenceStrings.find(seq.getAligned());
				
				if (itStrings == sequenceStrings.end()) { //this is a new unique sequence
					//output to unique fasta file
					seq.printSequence(outFasta);
					
					if (oldNameMapFName != "") {
						itNames = nameMap.find(seq.getName());
						
						if (itNames == nameMap.end()) { //namefile and fastafile do not match
							m->mothurOut("[ERROR]: " + seq.getName() + " is in your fasta file, and not in your namefile, please correct."); m->mothurOutEndLine();
						}else {
							sequenceStrings[seq.getAligned()] = itNames->second;
							nameFileOrder.push_back(seq.getAligned());
						}
					}else if (countfile != "") { 
                        ct.getNumSeqs(seq.getName()); //checks to make sure seq is in table
                        sequenceStrings[seq.getAligned()] = seq.getName();	nameFileOrder.push_back(seq.getAligned());
                    }else {	sequenceStrings[seq.getAligned()] = seq.getName();	nameFileOrder.push_back(seq.getAligned()); }
				}else { //this is a dup
					if (oldNameMapFName != "") {
						itNames = nameMap.find(seq.getName());
						
						if (itNames == nameMap.end()) { //namefile and fastafile do not match
							m->mothurOut("[ERROR]: " + seq.getName() + " is in your fasta file, and not in your namefile, please correct."); m->mothurOutEndLine();
						}else {
							sequenceStrings[seq.getAligned()] += "," + itNames->second;
						}
                    }else if (countfile != "") { 
                        int num = ct.getNumSeqs(seq.getName()); //checks to make sure seq is in table
                        if (num != 0) { //its in the table
                            ct.mergeCounts(itStrings->second, seq.getName()); //merges counts and saves in uniques name
                        }
                    }else {	sequenceStrings[seq.getAligned()] += "," + seq.getName();	}
				}
				
				count++;
			}
			
			m->gobble(in);
			
			if(count % 1000 == 0)	{ m->mothurOut(toString(count) + "\t" + toString(sequenceStrings.size())); m->mothurOutEndLine();	}
		}
		
		if(count % 1000 != 0)	{ m->mothurOut(toString(count) + "\t" + toString(sequenceStrings.size())); m->mothurOutEndLine();	}
		
		in.close();
		outFasta.close();
		
		if (m->control_pressed) { m->mothurRemove(outFastaFile); return 0; }
		
		//print new names file
		ofstream outNames;
		if (countfile == "") { m->openOutputFile(outNameFile, outNames); outputNames.push_back(outNameFile); outputTypes["name"].push_back(outNameFile);  }
        else { m->openOutputFile(outCountFile, outNames); ct.printHeaders(outNames); outputTypes["count"].push_back(outCountFile); outputNames.push_back(outCountFile); }
		
		for (int i = 0; i < nameFileOrder.size(); i++) {
			if (m->control_pressed) { outputTypes.clear(); m->mothurRemove(outFastaFile); outNames.close(); for (int j = 0; j < outputNames.size(); j++) { m->mothurRemove(outputNames[j]); } return 0; }
			
			itStrings = sequenceStrings.find(nameFileOrder[i]);
			
			if (itStrings != sequenceStrings.end()) {
                if (countfile == "") {
                    //get rep name
                    int pos = (itStrings->second).find_first_of(',');
                    
                    if (pos == string::npos) { // only reps itself
                        outNames << itStrings->second << '\t' << itStrings->second << endl;
                    }else {
                        outNames << (itStrings->second).substr(0, pos) << '\t' << itStrings->second << endl;
                    }
                }else {  ct.printSeq(outNames, itStrings->second);  }
			}else{ m->mothurOut("[ERROR]: mismatch in namefile print."); m->mothurOutEndLine(); m->control_pressed = true; }
		}
		outNames.close();
		
		if (m->control_pressed) { outputTypes.clear(); m->mothurRemove(outFastaFile); for (int j = 0; j < outputNames.size(); j++) { m->mothurRemove(outputNames[j]); }  return 0; }
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		outputNames.push_back(outFastaFile);   outputTypes["fasta"].push_back(outFastaFile);  
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
        
        itTypes = outputTypes.find("count");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setCountTableFile(current); }
		}
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "DeconvoluteCommand", "execute");
		exit(1);
	}
}
/**************************************************************************************/
