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
#include "fastqread.h"

//**********************************************************************************************************************
vector<string> MakeFastQCommand::setParameters(){	
	try {
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none","fastq",false,true,true); parameters.push_back(pfasta);
		CommandParameter pqfile("qfile", "InputTypes", "", "", "none", "none", "none","fastq",false,true,true); parameters.push_back(pqfile);
		CommandParameter pformat("format", "Multiple", "sanger-illumina-solexa-illumina1.8+", "illumina1.8+", "", "", "","",false,false,true); parameters.push_back(pformat);
        CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
        abort = false; calledHelp = false;
        
        //initialize outputTypes
        vector<string> tempOutNames;
        outputTypes["fastq"] = tempOutNames;
       
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
		helpString += "The format parameter is used to indicate whether your sequences are sanger, solexa, illumina1.8+ or illumina, default=illumina1.8+.\n";
		helpString += "The make.fastq command should be in the following format: make.fastq(qfile=yourQualityFile, fasta=yourFasta).\n";
		helpString += "Example make.fastq(fasta=amazon.fasta, qfile=amazon.qual).\n";
		;
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
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeFastQCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
MakeFastQCommand::MakeFastQCommand(string option) : Command()  {
	try {

		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }
		
		else {
			OptionParser parser(option, setParameters());
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			fastafile = validParameter.validFile(parameters, "fasta");
			if (fastafile == "not open") { abort = true; fastafile = ""; }
			else if (fastafile == "not found") {  		
				fastafile = current->getFastaFile(); 
				if (fastafile != "") {  m->mothurOut("Using " + fastafile + " as input file for the fasta parameter.\n");  }
				else { 	m->mothurOut("You have no current fastafile and the fasta parameter is required.\n");  abort = true; }
			}else { current->setFastaFile(fastafile); }	
			
			qualfile = validParameter.validFile(parameters, "qfile");
			if (qualfile == "not open") { abort = true; qualfile = ""; }
			else if (qualfile == "not found") {  			
				qualfile = current->getQualFile(); 
				if (qualfile != "") {  m->mothurOut("Using " + qualfile + " as input file for the qfile parameter.\n");  }
				else { 	m->mothurOut("You have no current qualfile and the qfile parameter is required.\n");  abort = true; }
			}else { current->setQualFile(qualfile); }	
			
			 
					if (outputdir == ""){    outputdir = util.hasPath(fastafile);		}
            
            format = validParameter.valid(parameters, "format");		if (format == "not found"){	format = "illumina1.8+";	}
            
            if ((format != "sanger") && (format != "illumina") && (format != "illumina1.8+") && (format != "solexa"))  {
                m->mothurOut(format + " is not a valid format. Your format choices are sanger, solexa, illumina1.8+ and illumina, aborting." ); m->mothurOutEndLine();
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
		
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
        map<string, string> variables; 
        variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(fastafile));
		string outputFile = getOutputFileName("fastq",variables);
		outputNames.push_back(outputFile); outputTypes["fastq"].push_back(outputFile);
		
		ofstream out;
		util.openOutputFile(outputFile, out);
		
		ifstream qFile;
		util.openInputFile(qualfile, qFile);
		
		ifstream fFile;
		util.openInputFile(fastafile, fFile);
        
		while (!fFile.eof() && !qFile.eof()) {
			
			if (m->getControl_pressed()) { break; }
			
			Sequence currSeq(fFile); gobble(fFile);
			QualityScores currQual(qFile);  gobble(qFile);
            
            FastqRead fread(currSeq, currQual, format);
            fread.printFastq(out);
		}
		
		fFile.close();
		qFile.close();
		out.close();
		
		if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); } return 0; }
		
		m->mothurOut("\nOutput File Names: \n"); 
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i] +"\n"); 	} m->mothurOutEndLine();
					
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




		
