/*
 *  makefastqcommand.cpp
 *  mothur
 *
 *  Created by westcott on 2/14/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "makefastqcommand.h"
#include "sequence.hpp"
#include "qualityscores.h"

//**********************************************************************************************************************
vector<string> MakeFastQCommand::setParameters(){	
	try {
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none","fastq",false,true,true); parameters.push_back(pfasta);
		CommandParameter pqfile("qfile", "InputTypes", "", "", "none", "none", "none","fastq",false,true,true); parameters.push_back(pqfile);
		CommandParameter pformat("format", "Multiple", "sanger-illumina-illumina1.8+", "sanger", "", "", "","",false,false); parameters.push_back(pformat);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "MakeFastQCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string MakeFastQCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The make.fastq command reads a fasta and quality file and creates a fastq file.\n";
		helpString += "The make.fastq command parameters are fasta, qfile and format.  fasta and qfile are required.\n";
		helpString += "The format parameter is used to indicate whether your sequences are sanger, illumina1.8+ or illumina, default=sanger.\n";
		helpString += "The make.fastq command should be in the following format: make.fastq(qfile=yourQualityFile, fasta=yourFasta).\n";
		helpString += "Example make.fastq(fasta=amazon.fasta, qfile=amazon.qual).\n";
		helpString += "Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFasta).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "MakeFastQCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string MakeFastQCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "fastq") {  pattern = "[filename],fastq"; } 
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeFastQCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
MakeFastQCommand::MakeFastQCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["fastq"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "MakeFastQCommand", "MakeFastQCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
MakeFastQCommand::MakeFastQCommand(string option)  {
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
			map<string,string>::iterator it;
			
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["fastq"] = tempOutNames;
			
						
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
				
				it = parameters.find("qfile");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["qfile"] = inputDir + it->second;		}
				}
			}
			
			
			//check for required parameters
			fastafile = validParameter.validFile(parameters, "fasta", true);
			if (fastafile == "not open") { abort = true; fastafile = ""; }
			else if (fastafile == "not found") {  		
				fastafile = m->getFastaFile(); 
				if (fastafile != "") {  m->mothurOut("Using " + fastafile + " as input file for the fasta parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current fastafile and the fasta parameter is required."); m->mothurOutEndLine(); abort = true; }
			}else { m->setFastaFile(fastafile); }	
			
			qualfile = validParameter.validFile(parameters, "qfile", true);
			if (qualfile == "not open") { abort = true; qualfile = ""; }
			else if (qualfile == "not found") {  			
				qualfile = m->getQualFile(); 
				if (qualfile != "") {  m->mothurOut("Using " + qualfile + " as input file for the qfile parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current qualfile and the qfile parameter is required."); m->mothurOutEndLine(); abort = true; }
			}else { m->setQualFile(qualfile); }	
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = m->hasPath(fastafile);		}
            
            format = validParameter.validFile(parameters, "format", false);		if (format == "not found"){	format = "sanger";	}
            
            if ((format != "sanger") && (format != "illumina") && (format != "illumina1.8+"))  { 
				m->mothurOut(format + " is not a valid format. Your format choices are sanger, illumina1.8+ and illumina, aborting." ); m->mothurOutEndLine();
				abort=true;
			}


		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "MakeFastQCommand", "MakeFastQCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int MakeFastQCommand::execute(){
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
        map<string, string> variables; 
        variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(fastafile));
		string outputFile = getOutputFileName("fastq",variables);
		outputNames.push_back(outputFile); outputTypes["fastq"].push_back(outputFile);
		
		ofstream out;
		m->openOutputFile(outputFile, out);
		
		ifstream qFile;
		m->openInputFile(qualfile, qFile);
		
		ifstream fFile;
		m->openInputFile(fastafile, fFile);
        
		while (!fFile.eof() && !qFile.eof()) {
			
			if (m->control_pressed) { break; }
			
			Sequence currSeq(fFile); m->gobble(fFile);
			QualityScores currQual(qFile);  m->gobble(qFile);
			
			if (currSeq.getName() != currQual.getName()) { m->mothurOut("[ERROR]: mismatch between fasta and quality files. Found " + currSeq.getName() + " in fasta file and " + currQual.getName() + " in quality file."); m->mothurOutEndLine(); m->control_pressed = true; }
			else {
				//print sequence
				out << '@' << currSeq.getName() << endl << currSeq.getAligned() << endl;
				
				string qualityString = convertQual(currQual.getQualityScores());
				
				//print quality info
				out << '+' << currQual.getName() << endl << qualityString << endl;
			}
			
		}
		
		fFile.close();
		qFile.close();
		out.close();
		
		if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); } return 0; }
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
					
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "MakeFastQCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
string MakeFastQCommand::convertQual(vector<int> qual) {
	try {
		string qualScores;
		
        for (int i = 0; i < qual.size(); i++) { 
            int controlChar = int('!');
            if (format == "illumina") {  controlChar = int('@');  }
            
			int temp = qual[i] + controlChar;
			char qualChar = (char) temp;
			
			qualScores += qualChar;
		}
		
		return qualScores;
	}
	catch(exception& e) {
		m->errorOut(e, "MakeFastQCommand", "convertQual");
		exit(1);
	}
}
//**********************************************************************************************************************




		
