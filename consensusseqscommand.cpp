/*
 *  consensusseqscommand.cpp
 *  Mothur
 *
 *  Created by westcott on 11/23/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "consensusseqscommand.h"
#include "sequence.hpp"
#include "inputdata.h"

//**********************************************************************************************************************
vector<string> ConsensusSeqsCommand::setParameters(){	
	try {
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(pfasta);
		CommandParameter pname("name", "InputTypes", "", "", "none", "none", "none",false,false); parameters.push_back(pname);
		CommandParameter plist("list", "InputTypes", "", "", "none", "none", "none",false,false); parameters.push_back(plist);
		CommandParameter plabel("label", "String", "", "", "", "", "",false,false); parameters.push_back(plabel);
		CommandParameter pcutoff("cutoff", "Number", "", "100", "", "", "",false,false); parameters.push_back(pcutoff);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ConsensusSeqsCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string ConsensusSeqsCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The consensus.seqs command can be used in 2 ways: create a consensus sequence from a fastafile, or with a listfile create a consensus sequence for each otu. Sequences must be aligned.\n";
		helpString += "The consensus.seqs command parameters are fasta, list, name, cutoff and label.\n";
		helpString += "The fasta parameter allows you to enter the fasta file containing your sequences, and is required, unless you have a valid current fasta file. \n";
		helpString += "The list parameter allows you to enter a your list file. \n";
		helpString += "The name parameter allows you to enter a names file associated with the fasta file. \n";
		helpString += "The label parameter allows you to select what distance levels you would like output files for, and are separated by dashes.\n";
		helpString += "The cutoff parameter allows you set a percentage of sequences that support the base. For example: cutoff=97 would only return a sequence that only showed ambiguities for bases that were not supported by at least 97% of sequences.\n";
		helpString += "The consensus.seqs command should be in the following format: \n";
		helpString += "consensus.seqs(fasta=yourFastaFile, list=yourListFile) \n";	
		helpString += "Example: consensus.seqs(fasta=abrecovery.align, list=abrecovery.fn.list) \n";
		helpString += "Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFastaFile).\n";	
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "ConsensusSeqsCommand", "getHelpString");
		exit(1);
	}
}

//**********************************************************************************************************************
ConsensusSeqsCommand::ConsensusSeqsCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["fasta"] = tempOutNames;
		outputTypes["name"] = tempOutNames;
		outputTypes["summary"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "ConsensusSeqsCommand", "ConsensusSeqsCommand");
		exit(1);
	}
}
//***************************************************************************************************************
ConsensusSeqsCommand::ConsensusSeqsCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		allLines = 1;
		
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
			outputTypes["fasta"] = tempOutNames;
			outputTypes["name"] = tempOutNames;
			outputTypes["summary"] = tempOutNames;
			
						
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
				
				it = parameters.find("list");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["list"] = inputDir + it->second;		}
				}
			}
			
			
			//check for parameters
			fastafile = validParameter.validFile(parameters, "fasta", true);
			if (fastafile == "not open") { abort = true; }
			else if (fastafile == "not found") { 			
				fastafile = m->getFastaFile(); 
				if (fastafile != "") { m->mothurOut("Using " + fastafile + " as input file for the fasta parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current fastafile and the fasta parameter is required."); m->mothurOutEndLine(); abort = true; }
			}	
			
			namefile = validParameter.validFile(parameters, "name", true);
			if (namefile == "not open") { abort = true; }
			else if (namefile == "not found") { namefile = ""; }	
			
			listfile = validParameter.validFile(parameters, "list", true);
			if (listfile == "not open") { abort = true; }
			else if (listfile == "not found") { listfile = "";  }	
			
			label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  m->splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
			
			string temp = validParameter.validFile(parameters, "cutoff", false);  if (temp == "not found") { temp = "100"; }
			convert(temp, cutoff); 
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = m->hasPath(fastafile);	}

		}
	}
	catch(exception& e) {
		m->errorOut(e, "ConsensusSeqsCommand", "ConsensusSeqsCommand");
		exit(1);
	}
}
//***************************************************************************************************************

int ConsensusSeqsCommand::execute(){
	try{
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		readFasta();
		
		if (m->control_pressed) { return 0; }
		
		if (namefile != "") { readNames(); }
		
		if (m->control_pressed) { return 0; }
		
				
		if (listfile == "") {
			
			ofstream outSummary;
			string outputSummaryFile = outputDir + m->getRootName(m->getSimpleName(fastafile)) + "cons.summary";
			m->openOutputFile(outputSummaryFile, outSummary);
			outSummary.setf(ios::fixed, ios::floatfield); outSummary.setf(ios::showpoint);
			outputNames.push_back(outputSummaryFile); outputTypes["summary"].push_back(outputSummaryFile);
			
			outSummary << "PositioninAlignment\tA\tT\tG\tC\tGap\tNumberofSeqs\tConsensusBase" << endl;
			
			ofstream outFasta;
			string outputFastaFile = outputDir + m->getRootName(m->getSimpleName(fastafile)) + "cons.fasta";
			m->openOutputFile(outputFastaFile, outFasta);
			outputNames.push_back(outputFastaFile); outputTypes["fasta"].push_back(outputFastaFile);
			
			vector<string> seqs;
			int seqLength = 0;
			for (map<string, string>::iterator it = nameMap.begin(); it != nameMap.end(); it++) {
				
				if (m->control_pressed) { outSummary.close(); outFasta.close(); for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); } return 0; }
				
				string seq = fastaMap[it->second];
				seqs.push_back(seq);
				
				if (seqLength == 0) { seqLength = seq.length(); }
				else if (seqLength != seq.length()) { m->mothurOut("[ERROR]: sequence are not the same length, please correct."); m->mothurOutEndLine(); m->control_pressed = true; }

			}
			
			vector< vector<float> > percentages; percentages.resize(5);
			for (int j = 0; j < percentages.size(); j++) { percentages[j].resize(seqLength, 0.0); }
			
			string consSeq = "";
			//get counts
			for (int j = 0; j < seqLength; j++) {
				
				if (m->control_pressed) { outSummary.close(); outFasta.close(); for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); } return 0; }
				
				vector<int> counts; counts.resize(5, 0); //A,T,G,C,Gap
				int numDots = 0;
				
				for (int i = 0; i < seqs.size(); i++) {
					
					if (seqs[i][j] == '.') { numDots++; }
					
					char base = toupper(seqs[i][j]);
					if (base == 'A') { counts[0]++; }
					else if (base == 'T') { counts[1]++; }
					else if (base == 'G') { counts[2]++; }
					else if (base == 'C') { counts[3]++; }
					else { counts[4]++; }
				}
				
				char conBase = '.';
				if (numDots != seqs.size()) { conBase = getBase(counts, seqs.size()); }
				
				consSeq += conBase;
				
				percentages[0][j] = counts[0] / (float) seqs.size();
				percentages[1][j] = counts[1] / (float) seqs.size();
				percentages[2][j] = counts[2] / (float) seqs.size();
				percentages[3][j] = counts[3] / (float) seqs.size();
				percentages[4][j] = counts[4] / (float) seqs.size();
				
			}
			
			for (int j = 0; j < seqLength; j++) { 
				outSummary << (j+1) << '\t' << percentages[0][j] << '\t'<< percentages[1][j] << '\t'<< percentages[2][j] << '\t' << percentages[3][j] << '\t' << percentages[4][j] << '\t' << seqs.size() << '\t' << consSeq[j] << endl;
			}
			
				
			outFasta << ">conseq" << endl << consSeq << endl;
			
			outSummary.close(); outFasta.close();
			
		
		}else {
			
						
			InputData* input = new InputData(listfile, "list");
			ListVector* list = input->getListVector();
			
			string lastLabel = list->getLabel();
			set<string> processedLabels;
			set<string> userLabels = labels;

			//as long as you are not at the end of the file or done wih the lines you want
			while((list != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
				
				if (m->control_pressed) {  for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); } delete list; delete input;  return 0;  }
				
				if(allLines == 1 || labels.count(list->getLabel()) == 1){			
					
					m->mothurOut(list->getLabel()); m->mothurOutEndLine();
					
					processList(list);
					
					processedLabels.insert(list->getLabel());
					userLabels.erase(list->getLabel());
				}
				
				if ((m->anyLabelsToProcess(list->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
					string saveLabel = list->getLabel();
					
					delete list; 
					
					list = input->getListVector(lastLabel);
					m->mothurOut(list->getLabel()); m->mothurOutEndLine();
					
					processList(list);
					
					processedLabels.insert(list->getLabel());
					userLabels.erase(list->getLabel());
					
					//restore real lastlabel to save below
					list->setLabel(saveLabel);
				}
				
				lastLabel = list->getLabel();
				
				delete list; list = NULL;
				
				//get next line to process
				list = input->getListVector();				
			}
			
			
			if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); } if (list != NULL) { delete list; } delete input; return 0;  }
			
			//output error messages about any remaining user labels
			set<string>::iterator it;
			bool needToRun = false;
			for (it = userLabels.begin(); it != userLabels.end(); it++) {  
				m->mothurOut("Your file does not include the label " + *it); 
				if (processedLabels.count(lastLabel) != 1) {
					m->mothurOut(". I will use " + lastLabel + "."); m->mothurOutEndLine();
					needToRun = true;
				}else {
					m->mothurOut(". Please refer to " + lastLabel + "."); m->mothurOutEndLine();
				}
			}
			
			//run last label if you need to
			if (needToRun == true)  {
				if (list != NULL) { delete list; }
				
				list = input->getListVector(lastLabel);
				
				m->mothurOut(list->getLabel()); m->mothurOutEndLine();
				
				processList(list);
				
				delete list; list = NULL;
			}
			
			if (list != NULL) { delete list; }
			delete input;
		}
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Name: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}	
		m->mothurOutEndLine();
		
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "ConsensusSeqsCommand", "execute");
		exit(1);
	}
}
//***************************************************************************************************************

int ConsensusSeqsCommand::processList(ListVector*& list){
	try{
		
		ofstream outSummary;
		string outputSummaryFile = outputDir + m->getRootName(m->getSimpleName(fastafile)) + list->getLabel() + ".cons.summary";
		m->openOutputFile(outputSummaryFile, outSummary);
		outSummary.setf(ios::fixed, ios::floatfield); outSummary.setf(ios::showpoint);
		outputNames.push_back(outputSummaryFile); outputTypes["summary"].push_back(outputSummaryFile);
		
		ofstream outName;
		string outputNameFile = outputDir + m->getRootName(m->getSimpleName(fastafile)) + list->getLabel() + ".cons.names";
		m->openOutputFile(outputNameFile, outName);
		outputNames.push_back(outputNameFile); outputTypes["name"].push_back(outputNameFile);
		
		ofstream outFasta;
		string outputFastaFile = outputDir + m->getRootName(m->getSimpleName(fastafile)) + list->getLabel() + ".cons.fasta";
		m->openOutputFile(outputFastaFile, outFasta);
		outputNames.push_back(outputFastaFile); outputTypes["fasta"].push_back(outputFastaFile);
		
		outSummary << "OTU#\tPositioninAlignment\tA\tT\tG\tC\tGap\tNumberofSeqs\tConsensusBase" << endl;
		
		for (int i = 0; i < list->getNumBins(); i++) {
			
			if (m->control_pressed) { outSummary.close(); outName.close(); outFasta.close(); return 0; }
			
			string bin = list->get(i);
			
			string newName = "";
			string consSeq = getConsSeq(bin, outSummary, newName, i);
			
			outFasta << ">seq" << (i+1) << endl << consSeq << endl;
			outName << "seq" << (i+1) << '\t' << "seq" << (i+1) << "," << newName << endl;
		}
		
		outSummary.close(); outName.close(); outFasta.close();
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "ConsensusSeqsCommand", "processList");
		exit(1);
	}
}

//***************************************************************************************************************
//made this smart enough to owrk with unique or non unique list file
string ConsensusSeqsCommand::getConsSeq(string bin, ofstream& outSummary, string& name, int binNumber){
	try{
		
		string consSeq = "";
		bool error = false;
		
		//the whole bin is the second column if no names file, otherwise build it
		name = bin;
		if (namefile != "") { name = ""; }
		
		vector<string> binNames;
		m->splitAtComma(bin, binNames);
		
		//get sequence strings for each name in the bin
		vector<string> seqs;
		
		set<string> addedAlready;
		int seqLength = 0;
		for (int i = 0; i < binNames.size(); i++) {
			
			map<string, string>::iterator it;
			
			it = nameMap.find(binNames[i]);
			if (it == nameMap.end()) { 
				if (namefile == "") { m->mothurOut("[ERROR]: " + binNames[i] + " is not in your fasta file, please correct."); m->mothurOutEndLine(); error = true; }
				else { m->mothurOut("[ERROR]: " + binNames[i] + " is not in your fasta or name file, please correct."); m->mothurOutEndLine(); error = true; }
				break;
			}else {
				
				//add sequence string to seqs vector to process below
				string seq = fastaMap[it->second];
				seqs.push_back(seq);
				
				if (seqLength == 0) { seqLength = seq.length(); }
				else if (seqLength != seq.length()) { m->mothurOut("[ERROR]: sequence are not the same length, please correct."); m->mothurOutEndLine(); error = true; break; }
				
				if (namefile != "") { 
					//did we add this line from name file already?
					if (addedAlready.count(it->second) == 0) {
						name += "," + nameFileMap[it->second];
						addedAlready.insert(it->second);
					}
				}
				
			}
		}
		
		if (error) { m->control_pressed = true; return consSeq; }
		
		if (namefile != "") { name = name.substr(1); }
		
		vector< vector<float> > percentages; percentages.resize(5);
		for (int j = 0; j < percentages.size(); j++) { percentages[j].resize(seqLength, 0.0); }
		
		//get counts
		for (int j = 0; j < seqLength; j++) {
			
			if (m->control_pressed) { return consSeq; }
			
			vector<int> counts; counts.resize(5, 0); //A,T,G,C,Gap
			int numDots = 0;
			
			for (int i = 0; i < seqs.size(); i++) {
				
				if (seqs[i][j] == '.') { numDots++; }
				
				char base = toupper(seqs[i][j]);
				if (base == 'A') { counts[0]++; }
				else if (base == 'T') { counts[1]++; }
				else if (base == 'G') { counts[2]++; }
				else if (base == 'C') { counts[3]++; }
				else { counts[4]++; }
			}
			
			char conBase = '.';
			if (numDots != seqs.size()) { conBase = getBase(counts, seqs.size()); }
			
			consSeq += conBase;
			
			percentages[0][j] = counts[0] / (float) seqs.size();
			percentages[1][j] = counts[1] / (float) seqs.size();
			percentages[2][j] = counts[2] / (float) seqs.size();
			percentages[3][j] = counts[3] / (float) seqs.size();
			percentages[4][j] = counts[4] / (float) seqs.size();
			
		}
		
		for (int j = 0; j < seqLength; j++) { 
			outSummary << (binNumber + 1) << '\t' << (j+1) << '\t' << percentages[0][j] << '\t'<< percentages[1][j] << '\t'<< percentages[2][j] << '\t' << percentages[3][j] << '\t' << percentages[4][j] << '\t' << seqs.size() << '\t' << consSeq[j] << endl;
		}
		
		return consSeq;
		
	}
	catch(exception& e) {
		m->errorOut(e, "ConsensusSeqsCommand", "getConsSeq");
		exit(1);
	}
}
//***************************************************************************************************************

char ConsensusSeqsCommand::getBase(vector<int> counts, int size){  //A,T,G,C,Gap
	try{
		/* A = adenine
		* C = cytosine
		* G = guanine
		* T = thymine
		* R = G A (purine)
		* Y = T C (pyrimidine)
		* K = G T (keto)
		* M = A C (amino)
		* S = G C (strong bonds)
		* W = A T (weak bonds)
		* B = G T C (all but A)
		* D = G A T (all but C)
		* H = A C T (all but G)
		* V = G C A (all but T)
		* N = A G C T (any) */
		
		char conBase = 'N';
		
		//zero out counts that don't make the cutoff
		float percentage = (100.0 - cutoff) / 100.0;
		int zeroCutoff = percentage * size;
		
		for (int i = 0; i < counts.size(); i++) {
			if (counts[i] < zeroCutoff) { counts[i] = 0; }
		}
		
		//any
		if ((counts[0] != 0) && (counts[1] != 0) && (counts[2] != 0) && (counts[3] != 0) && (counts[4] != 0)) {  conBase = 'n'; }
		//any no gap
		else if ((counts[0] != 0) && (counts[1] != 0) && (counts[2] != 0) && (counts[3] != 0) && (counts[4] == 0)) {  conBase = 'N'; }
		//all but T
		else if ((counts[0] != 0) && (counts[1] == 0) && (counts[2] != 0) && (counts[3] != 0) && (counts[4] != 0)) {  conBase = 'v'; }	
		//all but T no gap
		else if ((counts[0] != 0) && (counts[1] == 0) && (counts[2] != 0) && (counts[3] != 0) && (counts[4] == 0)) {  conBase = 'V'; }	
		//all but G
		else if ((counts[0] != 0) && (counts[1] != 0) && (counts[2] == 0) && (counts[3] != 0) && (counts[4] != 0)) {  conBase = 'h'; }	
		//all but G no gap
		else if ((counts[0] != 0) && (counts[1] != 0) && (counts[2] == 0) && (counts[3] != 0) && (counts[4] == 0)) {  conBase = 'H'; }	
		//all but C
		else if ((counts[0] != 0) && (counts[1] != 0) && (counts[2] != 0) && (counts[3] == 0) && (counts[4] != 0)) {  conBase = 'd'; }	
		//all but C no gap
		else if ((counts[0] != 0) && (counts[1] != 0) && (counts[2] != 0) && (counts[3] == 0) && (counts[4] == 0)) {  conBase = 'D'; }	
		//all but A
		else if ((counts[0] == 0) && (counts[1] != 0) && (counts[2] != 0) && (counts[3] != 0) && (counts[4] != 0)) {  conBase = 'b'; }	
		//all but A no gap
		else if ((counts[0] == 0) && (counts[1] != 0) && (counts[2] != 0) && (counts[3] != 0) && (counts[4] == 0)) {  conBase = 'B'; }	
		//W = A T (weak bonds)
		else if ((counts[0] != 0) && (counts[1] != 0) && (counts[2] == 0) && (counts[3] == 0) && (counts[4] != 0)) {  conBase = 'w'; }	
		//W = A T (weak bonds) no gap
		else if ((counts[0] != 0) && (counts[1] != 0) && (counts[2] == 0) && (counts[3] == 0) && (counts[4] == 0)) {  conBase = 'W'; }	
		//S = G C (strong bonds)
		else if ((counts[0] == 0) && (counts[1] == 0) && (counts[2] != 0) && (counts[3] != 0) && (counts[4] != 0)) {  conBase = 's'; }	
		//S = G C (strong bonds) no gap
		else if ((counts[0] == 0) && (counts[1] == 0) && (counts[2] != 0) && (counts[3] != 0) && (counts[4] == 0)) {  conBase = 'S'; }	
		//M = A C (amino)
		else if ((counts[0] != 0) && (counts[1] == 0) && (counts[2] == 0) && (counts[3] != 0) && (counts[4] != 0)) {  conBase = 'm'; }	
		//M = A C (amino) no gap
		else if ((counts[0] != 0) && (counts[1] == 0) && (counts[2] == 0) && (counts[3] != 0) && (counts[4] == 0)) {  conBase = 'M'; }	
		//K = G T (keto)
		else if ((counts[0] == 0) && (counts[1] != 0) && (counts[2] != 0) && (counts[3] == 0) && (counts[4] != 0)) {  conBase = 'k'; }	
		//K = G T (keto) no gap
		else if ((counts[0] == 0) && (counts[1] != 0) && (counts[2] != 0) && (counts[3] == 0) && (counts[4] == 0)) {  conBase = 'K'; }	
		//Y = T C (pyrimidine)
		else if ((counts[0] == 0) && (counts[1] != 0) && (counts[2] == 0) && (counts[3] != 0) && (counts[4] != 0)) {  conBase = 'y'; }	
		//Y = T C (pyrimidine) no gap
		else if ((counts[0] == 0) && (counts[1] != 0) && (counts[2] == 0) && (counts[3] != 0) && (counts[4] == 0)) {  conBase = 'Y'; }	
		//R = G A (purine)
		else if ((counts[0] != 0) && (counts[1] == 0) && (counts[2] != 0) && (counts[3] == 0) && (counts[4] != 0)) {  conBase = 'r'; }	
		//R = G A (purine) no gap
		else if ((counts[0] != 0) && (counts[1] == 0) && (counts[2] != 0) && (counts[3] == 0) && (counts[4] == 0)) {  conBase = 'R'; }	
		//only A
		else if ((counts[0] != 0) && (counts[1] == 0) && (counts[2] == 0) && (counts[3] == 0) && (counts[4] != 0)) {  conBase = 'a'; }	
		//only A no gap
		else if ((counts[0] != 0) && (counts[1] == 0) && (counts[2] == 0) && (counts[3] == 0) && (counts[4] == 0)) {  conBase = 'A'; }	
		//only T
		else if ((counts[0] == 0) && (counts[1] != 0) && (counts[2] == 0) && (counts[3] == 0) && (counts[4] != 0)) {  conBase = 't'; }	
		//only T no gap
		else if ((counts[0] == 0) && (counts[1] != 0) && (counts[2] == 0) && (counts[3] == 0) && (counts[4] == 0)) {  conBase = 'T'; }	
		//only G
		else if ((counts[0] == 0) && (counts[1] == 0) && (counts[2] != 0) && (counts[3] == 0) && (counts[4] != 0)) {  conBase = 'g'; }	
		//only G no gap
		else if ((counts[0] == 0) && (counts[1] == 0) && (counts[2] != 0) && (counts[3] == 0) && (counts[4] == 0)) {  conBase = 'G'; }	
		//only C
		else if ((counts[0] == 0) && (counts[1] == 0) && (counts[2] == 0) && (counts[3] != 0) && (counts[4] != 0)) {  conBase = 'c'; }	
		//only C no gap
		else if ((counts[0] == 0) && (counts[1] == 0) && (counts[2] == 0) && (counts[3] != 0) && (counts[4] == 0)) {  conBase = 'C'; }	
		//only gap
		else if ((counts[0] == 0) && (counts[1] == 0) && (counts[2] == 0) && (counts[3] == 0) && (counts[4] != 0)) {  conBase = '-'; }
		//cutoff removed all counts
		else if ((counts[0] == 0) && (counts[1] == 0) && (counts[2] == 0) && (counts[3] == 0) && (counts[4] == 0)) {  conBase = 'N'; }
		else{ m->mothurOut("[ERROR]: cannot find consensus base."); m->mothurOutEndLine(); }
		
		return conBase;
		
	}
	catch(exception& e) {
		m->errorOut(e, "ConsensusSeqsCommand", "getBase");
		exit(1);
	}
}

//***************************************************************************************************************

int ConsensusSeqsCommand::readFasta(){
	try{
		
		ifstream in;
		m->openInputFile(fastafile, in);
		
		while (!in.eof()) {
			
			if (m->control_pressed) { break; }
			
			Sequence seq(in); m->gobble(in);
			string name = seq.getName();
			
			if (name != "") {
				fastaMap[name] = seq.getAligned();
				nameMap[name] = name; //set nameMap incase no names file
				nameFileMap[name] = name;
			}
		}
		
		in.close();
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "ConsensusSeqsCommand", "readFasta");
		exit(1);
	}
}
//***************************************************************************************************************

int ConsensusSeqsCommand::readNames(){
	 try{
		 
		 ifstream in;
		 m->openInputFile(namefile, in);
		 
		 string thisname, repnames;
		 map<string, string>::iterator it;
		 
		 bool error = false;
		 
		 while(!in.eof()){
			 
			 if (m->control_pressed) { break; }
			 
			 in >> thisname;		m->gobble(in);		//read from first column
			 in >> repnames;			//read from second column
			 
			 it = nameMap.find(thisname);
			 if (it != nameMap.end()) { //then this sequence was in the fastafile
				 
				 vector<string> splitRepNames;
				 m->splitAtComma(repnames, splitRepNames);
				 
				 nameFileMap[thisname] = repnames;	//for later when outputting the new namesFile if the list file is unique
				 for (int i = 0; i < splitRepNames.size(); i++) { nameMap[splitRepNames[i]] = thisname; }
				 
			 }else{	m->mothurOut("[ERROR]: " + thisname + " is not in the fasta file, please correct."); m->mothurOutEndLine(); error = true; }
			 
			 m->gobble(in);
		 }
		 
		 in.close();
		 
		 if (error) { m->control_pressed = true; }
 
		 return 0;
 
	}
	 catch(exception& e) {
		 m->errorOut(e, "ConsensusSeqsCommand", "readNames");
		 exit(1);
	 }
 }
 
//***************************************************************************************************************


