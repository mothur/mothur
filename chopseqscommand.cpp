/*
 *  chopseqscommand.cpp
 *  Mothur
 *
 *  Created by westcott on 5/10/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "chopseqscommand.h"
#include "sequence.hpp"

//**********************************************************************************************************************

ChopSeqsCommand::ChopSeqsCommand(string option)  {
	try {
		abort = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"fasta","end","outputdir","inputdir"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			map<string,string>::iterator it;
			
			//check to make sure all parameters are valid for command
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";		}
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("fasta");
				//user has given a template file
				if(it != parameters.end()){ 
					path = hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["fasta"] = inputDir + it->second;		}
				}
			}

			//check for required parameters
			fastafile = validParameter.validFile(parameters, "fasta", true);
			if (fastafile == "not open") { abort = true; }
			else if (fastafile == "not found") {  m->mothurOut("You must provide a fasta file."); m->mothurOutEndLine(); abort = true; }  	
			
			string temp = validParameter.validFile(parameters, "end", false);	
			if (temp == "not found") {  m->mothurOut("You must provide an end for the chops.seqs command."); m->mothurOutEndLine(); abort = true; } 
			else {	
				convert(temp, end);   
				if (end < 0) { m->mothurOut("End must be positive."); m->mothurOutEndLine(); abort = true;  }
			}
			
		}

	}
	catch(exception& e) {
		m->errorOut(e, "ChopSeqsCommand", "ChopSeqsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

void ChopSeqsCommand::help(){
	try {
		m->mothurOut("The chop.seqs command reads a fasta file and outputs a .chop.fasta with sequences trimmed to the end position.\n");
		m->mothurOut("The chop.seqs command parameters are fasta and end, both are required.\n");
		m->mothurOut("The chop.seqs command should be in the following format: chop.seqs(fasta=yourFasta, end=yourEnd).\n");
		m->mothurOut("Example chop.seqs(fasta=amazon.fasta, end=200).\n");
		m->mothurOut("Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFasta).\n\n");
	}
	catch(exception& e) {
		m->errorOut(e, "ChopSeqsCommand", "help");
		exit(1);
	}
}

//**********************************************************************************************************************

int ChopSeqsCommand::execute(){
	try {
		
		if (abort == true) { return 0; }
		
		string outputFileName = outputDir + getRootName(getSimpleName(fastafile)) + "chop.fasta";
		
		ofstream out;
		openOutputFile(outputFileName, out);
		
		ifstream in;
		openInputFile(fastafile, in);
		
		while (!in.eof()) {
		
			Sequence seq(in, "no align");
			
			if (m->control_pressed) { in.close(); out.close(); remove(outputFileName.c_str()); return 0;  }
			
			if (seq.getName() != "") {
				string temp = seq.getUnaligned();
				
				//output sequence name
				out << ">" << seq.getName() << endl;
				
				//if needed trim sequence
				if (temp.length() > end) {  temp = temp.substr(0, end);		}
				
				//output trimmed sequence	
				out << temp << endl;
			}
		}
		in.close();
		out.close();
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Name: "); m->mothurOutEndLine();
		m->mothurOut(outputFileName); m->mothurOutEndLine();	
		m->mothurOutEndLine();
		
		return 0;		
	}

	catch(exception& e) {
		m->errorOut(e, "ChopSeqsCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************


