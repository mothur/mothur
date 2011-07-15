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
vector<string> ChopSeqsCommand::setParameters(){	
	try {
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(pfasta);
		CommandParameter pnumbases("numbases", "Number", "", "0", "", "", "",false,true); parameters.push_back(pnumbases);
		CommandParameter pcountgaps("countgaps", "Boolean", "", "F", "", "", "",false,false); parameters.push_back(pcountgaps);
		CommandParameter pshort("short", "Boolean", "", "F", "", "", "",false,false); parameters.push_back(pshort);
		CommandParameter pkeep("keep", "Multiple", "front-back", "front", "", "", "",false,false); parameters.push_back(pkeep);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ChopSeqsCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string ChopSeqsCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The chop.seqs command reads a fasta file and outputs a .chop.fasta containing the trimmed sequences. Note: If a sequence is completely 'chopped', an accnos file will be created with the names of the sequences removed. \n";
		helpString += "The chop.seqs command parameters are fasta, numbases, countgaps and keep. fasta is required unless you have a valid current fasta file. numbases is required.\n";
		helpString += "The chop.seqs command should be in the following format: chop.seqs(fasta=yourFasta, numbases=yourNum, keep=yourKeep).\n";
		helpString += "The numbases parameter allows you to specify the number of bases you want to keep.\n";
		helpString += "The keep parameter allows you to specify whether you want to keep the front or the back of your sequence, default=front.\n";
		helpString += "The countgaps parameter allows you to specify whether you want to count gaps as bases, default=false.\n";
		helpString += "The short parameter allows you to specify you want to keep sequences that are too short to chop, default=false.\n";
		helpString += "For example, if you ran chop.seqs with numbases=200 and short=t, if a sequence had 100 bases mothur would keep the sequence rather than eliminate it.\n";
		helpString += "Example chop.seqs(fasta=amazon.fasta, numbases=200, keep=front).\n";
		helpString += "Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFasta).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "ChopSeqsCommand", "getHelpString");
		exit(1);
	}
}

//**********************************************************************************************************************
ChopSeqsCommand::ChopSeqsCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["fasta"] = tempOutNames;
		outputTypes["accnos"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "ChopSeqsCommand", "ChopSeqsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
ChopSeqsCommand::ChopSeqsCommand(string option)  {
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
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["fasta"] = tempOutNames;
			outputTypes["accnos"] = tempOutNames;
		
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
			}

			//check for required parameters
			fastafile = validParameter.validFile(parameters, "fasta", true);
			if (fastafile == "not open") { abort = true; }
			else if (fastafile == "not found") {  				//if there is a current fasta file, use it
				fastafile = m->getFastaFile(); 
				if (fastafile != "") { m->mothurOut("Using " + fastafile + " as input file for the fasta parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current fastafile and the fasta parameter is required."); m->mothurOutEndLine(); abort = true; }
			}else { m->setFastaFile(fastafile); } 	
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = m->hasPath(fastafile);	}
			
			string temp = validParameter.validFile(parameters, "numbases", false);	if (temp == "not found") { temp = "0"; } 
			convert(temp, numbases);   
			
			temp = validParameter.validFile(parameters, "countgaps", false);	if (temp == "not found") { temp = "f"; } 
			countGaps = m->isTrue(temp);  
			
			temp = validParameter.validFile(parameters, "short", false);	if (temp == "not found") { temp = "f"; } 
			Short = m->isTrue(temp);   
		
			keep = validParameter.validFile(parameters, "keep", false);		if (keep == "not found") { keep = "front"; } 
				
			if (numbases == 0)  { m->mothurOut("You must provide the number of bases you want to keep for the chops.seqs command."); m->mothurOutEndLine(); abort = true;  }
		}

	}
	catch(exception& e) {
		m->errorOut(e, "ChopSeqsCommand", "ChopSeqsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int ChopSeqsCommand::execute(){
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		string outputFileName = outputDir + m->getRootName(m->getSimpleName(fastafile)) + "chop.fasta";
		string outputFileNameAccnos = outputDir + m->getRootName(m->getSimpleName(fastafile)) + "chop.accnos";
		
		ofstream out;
		m->openOutputFile(outputFileName, out);
		
		ofstream outAcc;
		m->openOutputFile(outputFileNameAccnos, outAcc);
		
		ifstream in;
		m->openInputFile(fastafile, in);
		
		bool wroteAccnos = false;
		
		while (!in.eof()) {
			
			Sequence seq(in);
			
			if (m->control_pressed) { outputTypes.clear(); in.close(); out.close(); outAcc.close(); m->mothurRemove(outputFileName); m->mothurRemove(outputFileNameAccnos); return 0;  }
			
			if (seq.getName() != "") {
				string newSeqString = getChopped(seq);
				
				//output trimmed sequence
				if (newSeqString != "") {
					out << ">" << seq.getName() << endl << newSeqString << endl;
				}else{
					outAcc << seq.getName() << endl;
					wroteAccnos = true;
				}
			}
		}
		in.close();
		out.close();
		outAcc.close();
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Name: "); m->mothurOutEndLine();
		m->mothurOut(outputFileName); m->mothurOutEndLine();	outputNames.push_back(outputFileName); outputTypes["fasta"].push_back(outputFileName);
		
		if (wroteAccnos) { m->mothurOut(outputFileNameAccnos); m->mothurOutEndLine(); outputNames.push_back(outputFileNameAccnos); outputTypes["accnos"].push_back(outputFileNameAccnos); }
		else {  m->mothurRemove(outputFileNameAccnos);  }
		
		m->mothurOutEndLine();
		
		//set fasta file as new current fastafile
		string current = "";
		itTypes = outputTypes.find("fasta");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setFastaFile(current); }
		}
		
		if (wroteAccnos) { //set accnos file as new current accnosfile
			itTypes = outputTypes.find("accnos");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setAccnosFile(current); }
			}
		}
		
		
		return 0;		
	}

	catch(exception& e) {
		m->errorOut(e, "ChopSeqsCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
string ChopSeqsCommand::getChopped(Sequence seq) {
	try {
		string temp = seq.getAligned();
		string tempUnaligned = seq.getUnaligned();
		
		if (countGaps) {
			//if needed trim sequence
			if (keep == "front") {//you want to keep the beginning
				int tempLength = temp.length();

				if (tempLength > numbases) { //you have enough bases to remove some
				
					int stopSpot = 0;
					int numBasesCounted = 0;
					
					for (int i = 0; i < temp.length(); i++) {
						//eliminate N's
						if (toupper(temp[i]) == 'N') { temp[i] == '.'; }
						
						numBasesCounted++; 
						
						if (numBasesCounted >= numbases) { stopSpot = i; break; }
					}
					
					if (stopSpot == 0) { temp = ""; }
					else {  temp = temp.substr(0, stopSpot);  }
							
				}else { 
					if (!Short) { temp = ""; } //sequence too short
				}
			}else { //you are keeping the back
				int tempLength = temp.length();
				if (tempLength > numbases) { //you have enough bases to remove some
					
					int stopSpot = 0;
					int numBasesCounted = 0;
					
					for (int i = (temp.length()-1); i >= 0; i--) {
						//eliminate N's
						if (toupper(temp[i]) == 'N') { temp[i] == '.'; }
						
						numBasesCounted++; 

						if (numBasesCounted >= numbases) { stopSpot = i; break; }
					}
				
					if (stopSpot == 0) { temp = ""; }
					else {  temp = temp.substr(stopSpot+1);  }
				}else { 
					if (!Short) { temp = ""; } //sequence too short
				}
			}

		}else{
				
			//if needed trim sequence
			if (keep == "front") {//you want to keep the beginning
				int tempLength = tempUnaligned.length();

				if (tempLength > numbases) { //you have enough bases to remove some
					
					int stopSpot = 0;
					int numBasesCounted = 0;
					
					for (int i = 0; i < temp.length(); i++) {
						//eliminate N's
						if (toupper(temp[i]) == 'N') { 
							temp[i] == '.'; 
							tempLength--;
							if (tempLength < numbases) { stopSpot = 0; break; }
						}
						
						if(isalpha(temp[i])) { numBasesCounted++; }
						
						if (numBasesCounted >= numbases) { stopSpot = i; break; }
					}
					
					if (stopSpot == 0) { temp = ""; }
					else {  temp = temp.substr(0, stopSpot);  }
							
				}else { 
					if (!Short) { temp = ""; } //sequence too short
				}				
			}else { //you are keeping the back
				int tempLength = tempUnaligned.length();
				if (tempLength > numbases) { //you have enough bases to remove some
					
					int stopSpot = 0;
					int numBasesCounted = 0;
					
					for (int i = (temp.length()-1); i >= 0; i--) {
						//eliminate N's
						if (toupper(temp[i]) == 'N') { 
							temp[i] == '.'; 
							tempLength--;
							if (tempLength < numbases) { stopSpot = 0; break; }
						}
						
						if(isalpha(temp[i])) { numBasesCounted++; }

						if (numBasesCounted >= numbases) { stopSpot = i; break; }
					}
				
					if (stopSpot == 0) { temp = ""; }
					else {  temp = temp.substr(stopSpot+1);  }
				}else { 
					if (!Short) { temp = ""; } //sequence too short
				}
			}
		}
		
		return temp;
	}
	catch(exception& e) {
		m->errorOut(e, "ChopSeqsCommand", "getChopped");
		exit(1);
	}
}
//**********************************************************************************************************************


