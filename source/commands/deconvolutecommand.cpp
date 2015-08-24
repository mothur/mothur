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
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none","fasta-name",false,true,true); parameters.push_back(pfasta);
		CommandParameter pname("name", "InputTypes", "", "", "namecount", "none", "none","name",false,false,true); parameters.push_back(pname);
        CommandParameter pcount("count", "InputTypes", "", "", "namecount", "none", "none","count",false,false,true); parameters.push_back(pcount);
        CommandParameter pformat("format", "Multiple", "count-name", "name", "", "", "","",false,false, true); parameters.push_back(pformat);
        CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
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
		helpString += "The unique.seqs command parameters are fasta, name, count and format.  fasta is required, unless there is a valid current fasta file.\n";
        helpString += "The name parameter is used to provide an existing name file associated with the fasta file. \n";
        helpString += "The count parameter is used to provide an existing count file associated with the fasta file. \n";
        helpString += "The format parameter is used to indicate what type of file you want outputted.  Choices are name and count, default=name unless count file used then default=count.\n";
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
string DeconvoluteCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "fasta") {  pattern = "[filename],unique,[extension]"; } 
        else if (type == "name") {  pattern = "[filename],names-[filename],[tag],names"; } 
        else if (type == "count") {  pattern = "[filename],count_table-[filename],[tag],count_table"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "DeconvoluteCommand", "getOutputPattern");
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
			fastafile = validParameter.validFile(parameters, "fasta", true);
			if (fastafile == "not open") { abort = true; }
			else if (fastafile == "not found") {
				fastafile = m->getFastaFile();
				if (fastafile != "") { m->mothurOut("Using " + fastafile + " as input file for the fasta parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current fastafile and the fasta parameter is required."); m->mothurOutEndLine(); abort = true; }
			}else { m->setFastaFile(fastafile); }
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += m->hasPath(fastafile); //if user entered a file with a path then preserve it
			}
			
			namefile = validParameter.validFile(parameters, "name", true);
			if (namefile == "not open") { namefile = ""; abort = true; }
			else if (namefile == "not found"){	namefile = "";	}
			else { m->setNameFile(namefile); }
            
            countfile = validParameter.validFile(parameters, "count", true);
			if (countfile == "not open") { abort = true; countfile = ""; }	
			else if (countfile == "not found") { countfile = ""; }
			else { m->setCountTableFile(countfile); }
			
            if ((countfile != "") && (namefile != "")) { m->mothurOut("When executing a unique.seqs command you must enter ONLY ONE of the following: count or name."); m->mothurOutEndLine(); abort = true; }
			
            format = validParameter.validFile(parameters, "format", false);
            if(format == "not found"){
                if (countfile != "") { format = "count";    }
                else { format = "name";                     }
            }
            
            if ((format != "name") && (format != "count")) {
                m->mothurOut(format + " is not a valid format option. Options are count or name.");
                if (countfile == "") { m->mothurOut("I will use name.\n"); format = "name"; }
                else {  m->mothurOut("I will use count.\n"); format = "count"; }
            }
            
			if (countfile == "") {
                if (namefile == "") {
                    vector<string> files; files.push_back(fastafile);
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
        map<string, string> variables; 
        variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(fastafile));
        string outNameFile = getOutputFileName("name", variables);
        string outCountFile = getOutputFileName("count", variables);
        variables["[extension]"] = m->getExtension(fastafile);
		string outFastaFile = getOutputFileName("fasta", variables);
		
		map<string, string> nameMap;
		map<string, string>::iterator itNames;
		if (namefile != "")  {
            m->readNames(namefile, nameMap);
            if (namefile == outNameFile){
                //prepare filenames and open files
                map<string, string> mvariables;
                mvariables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(fastafile));
                mvariables["[tag]"] = "unique";
                outNameFile = getOutputFileName("name", mvariables);
            }
        }
        CountTable ct;
        if (countfile != "")  {  
            ct.readTable(countfile, true, false);
            if (countfile == outCountFile){
                //prepare filenames and open files
                map<string, string> mvariables;
                mvariables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(fastafile));
                mvariables["[tag]"] = "unique";
                outCountFile = getOutputFileName("count", mvariables);   }
        }
		
		if (m->control_pressed) { return 0; }
		
		ifstream in; 
		m->openInputFile(fastafile, in);
		
		ofstream outFasta;
		m->openOutputFile(outFastaFile, outFasta);
		
		map<string, string> sequenceStrings; //sequenceString -> list of names.  "atgc...." -> seq1,seq2,seq3.
		map<string, string>::iterator itStrings;
		set<string> nameInFastaFile; //for sanity checking
		set<string>::iterator itname;
		vector<string> nameFileOrder;
        CountTable newCt;
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
					
					if (namefile != "") {
						itNames = nameMap.find(seq.getName());
						
						if (itNames == nameMap.end()) { //namefile and fastafile do not match
							m->mothurOut("[ERROR]: " + seq.getName() + " is in your fasta file, and not in your namefile, please correct."); m->mothurOutEndLine();
						}else {
                            if (format == "name") { sequenceStrings[seq.getAligned()] = itNames->second;  nameFileOrder.push_back(seq.getAligned());
                            }else {  newCt.push_back(seq.getName(), m->getNumNames(itNames->second)); sequenceStrings[seq.getAligned()] = seq.getName();	nameFileOrder.push_back(seq.getAligned()); }
						}
					}else if (countfile != "") {
                        if (format == "name") {
                            int numSeqs = ct.getNumSeqs(seq.getName());
                            string expandedName = seq.getName()+"_0";
                            for (int i = 1; i < numSeqs; i++) {  expandedName += "," + seq.getName() + "_" + toString(i);  }
                            sequenceStrings[seq.getAligned()] = expandedName;  nameFileOrder.push_back(seq.getAligned());
                        }else {
                            ct.getNumSeqs(seq.getName()); //checks to make sure seq is in table
                            sequenceStrings[seq.getAligned()] = seq.getName();	nameFileOrder.push_back(seq.getAligned());
                        }
                    }else {
                        if (format == "name") { sequenceStrings[seq.getAligned()] = seq.getName();	nameFileOrder.push_back(seq.getAligned()); }
                        else {  newCt.push_back(seq.getName()); sequenceStrings[seq.getAligned()] = seq.getName();	nameFileOrder.push_back(seq.getAligned()); }
                    }
				}else { //this is a dup
					if (namefile != "") {
						itNames = nameMap.find(seq.getName());
						
						if (itNames == nameMap.end()) { //namefile and fastafile do not match
							m->mothurOut("[ERROR]: " + seq.getName() + " is in your fasta file, and not in your namefile, please correct."); m->mothurOutEndLine();
						}else {
                            if (format == "name") { sequenceStrings[seq.getAligned()] += "," + itNames->second;  }
                            else {  int currentReps = newCt.getNumSeqs(itStrings->second);  newCt.setNumSeqs(itStrings->second, currentReps+(m->getNumNames(itNames->second)));  }
						}
                    }else if (countfile != "") {
                        if (format == "name") {
                            int numSeqs = ct.getNumSeqs(seq.getName());
                            string expandedName = seq.getName()+"_0";
                            for (int i = 1; i < numSeqs; i++) {  expandedName += "," + seq.getName() + "_" + toString(i);  }
                            sequenceStrings[seq.getAligned()] += "," + expandedName;
                        }else {
                            int num = ct.getNumSeqs(seq.getName()); //checks to make sure seq is in table
                            if (num != 0) { //its in the table
                                ct.mergeCounts(itStrings->second, seq.getName()); //merges counts and saves in uniques name
                            }
                        }
                    }else {
                        if (format == "name") {  sequenceStrings[seq.getAligned()] += "," + seq.getName(); }
                        else {  int currentReps = newCt.getNumSeqs(itStrings->second); newCt.setNumSeqs(itStrings->second, currentReps+1);  }
                    }
				}
				count++;
			}
			
			m->gobble(in);
			
			if(count % 1000 == 0)	{ m->mothurOutJustToScreen(toString(count) + "\t" + toString(sequenceStrings.size()) + "\n");	}
		}
		
		if(count % 1000 != 0)	{ m->mothurOut(toString(count) + "\t" + toString(sequenceStrings.size())); m->mothurOutEndLine();	}
		
		in.close();
		outFasta.close();
		
		if (m->control_pressed) { m->mothurRemove(outFastaFile); return 0; }
        
		//print new names file
		ofstream outNames;
		if (format == "name") { m->openOutputFile(outNameFile, outNames); outputNames.push_back(outNameFile); outputTypes["name"].push_back(outNameFile);   }
        else { m->openOutputFile(outCountFile, outNames); outputTypes["count"].push_back(outCountFile); outputNames.push_back(outCountFile);                }
        
        if ((countfile != "") && (format == "count")) { ct.printHeaders(outNames); }
        else if ((countfile == "") && (format == "count")) { newCt.printHeaders(outNames); }
		
		for (int i = 0; i < nameFileOrder.size(); i++) {
			if (m->control_pressed) { outputTypes.clear(); m->mothurRemove(outFastaFile); outNames.close(); for (int j = 0; j < outputNames.size(); j++) { m->mothurRemove(outputNames[j]); } return 0; }
			
			itStrings = sequenceStrings.find(nameFileOrder[i]);
			
			if (itStrings != sequenceStrings.end()) {
                if (format == "name") {
                    //get rep name
                    int pos = (itStrings->second).find_first_of(',');
                    
                    if (pos == string::npos) { // only reps itself
                        outNames << itStrings->second << '\t' << itStrings->second << endl;
                    }else {
                        outNames << (itStrings->second).substr(0, pos) << '\t' << itStrings->second << endl;
                    }
                }else {
                    if (countfile != "") {  ct.printSeq(outNames, itStrings->second);  }
                    else if (format == "count")  {  newCt.printSeq(outNames, itStrings->second);  }
                }
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
