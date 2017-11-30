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
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none","fasta-name",false,true,true); parameters.push_back(pfasta);
		CommandParameter pname("name", "InputTypes", "", "", "namecount", "none", "none","name",false,false,true); parameters.push_back(pname);
        CommandParameter pcount("count", "InputTypes", "", "", "namecount", "none", "none","count",false,false,true); parameters.push_back(pcount);
		CommandParameter plist("list", "InputTypes", "", "", "none", "none", "none","fasta-name",false,false,true); parameters.push_back(plist);
		CommandParameter plabel("label", "String", "", "", "", "", "","",false,false); parameters.push_back(plabel);
		CommandParameter pcutoff("cutoff", "Number", "", "100", "", "", "","",false,false); parameters.push_back(pcutoff);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
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
		helpString += "The consensus.seqs command parameters are fasta, list, name, count, cutoff and label.\n";
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
string ConsensusSeqsCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "fasta") {  pattern = "[filename],cons.fasta-[filename],[tag],cons.fasta"; } 
        else if (type == "name") {  pattern = "[filename],cons.names-[filename],[tag],cons.names"; } 
        else if (type == "count") {  pattern = "[filename],cons.count_table-[filename],[tag],cons.count_table"; }
        else if (type == "summary") {  pattern = "[filename],cons.summary-[filename],[tag],cons.summary"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "ConsensusSeqsCommand", "getOutputPattern");
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
        outputTypes["count"] = tempOutNames;
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
            outputTypes["count"] = tempOutNames;
			outputTypes["summary"] = tempOutNames;
			
						
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.valid(parameters, "inputdir");	
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("fasta");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["fasta"] = inputDir + it->second;		}
				}
				
				it = parameters.find("name");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["name"] = inputDir + it->second;		}
				}
				
				it = parameters.find("list");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["list"] = inputDir + it->second;		}
				}
                
                it = parameters.find("count");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["count"] = inputDir + it->second;		}
				}
			}
			
			
			//check for parameters
			fastafile = validParameter.validFile(parameters, "fasta");
			if (fastafile == "not open") { abort = true; }
			else if (fastafile == "not found") { 			
				fastafile = current->getFastaFile(); 
				if (fastafile != "") { m->mothurOut("Using " + fastafile + " as input file for the fasta parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current fastafile and the fasta parameter is required."); m->mothurOutEndLine(); abort = true; }
			}else { current->setFastaFile(fastafile); }	
			
			namefile = validParameter.validFile(parameters, "name");
			if (namefile == "not open") { namefile = ""; abort = true; }
			else if (namefile == "not found") { namefile = ""; }
			else { current->setNameFile(namefile); }
			
            countfile = validParameter.validFile(parameters, "count");
			if (countfile == "not open") { abort = true; countfile = ""; }	
			else if (countfile == "not found") { countfile = ""; }
			else { current->setCountFile(countfile); }
			
            if ((countfile != "") && (namefile != "")) { m->mothurOut("You must enter ONLY ONE of the following: count or name."); m->mothurOutEndLine(); abort = true; }
            
			listfile = validParameter.validFile(parameters, "list");
			if (listfile == "not open") { abort = true; }
			else if (listfile == "not found") { listfile = "";  }	
			else { current->setListFile(listfile); }
			
			label = validParameter.valid(parameters, "label");			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  util.splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
			
			string temp = validParameter.valid(parameters, "cutoff");  if (temp == "not found") { temp = "100"; }
			util.mothurConvert(temp, cutoff); 
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){	outputDir = util.hasPath(fastafile);	}
			
            if (countfile == "") {
                if (namefile == ""){
                    vector<string> files; files.push_back(fastafile); 
                    if (!current->getMothurCalling())  {  parser.getNameFile(files);  }
                }
            }
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
		
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
        long start = time(NULL);
        
		readFasta();
		
		if (m->getControl_pressed()) { return 0; }
		
		if (namefile != "") { readNames(); }
        if (countfile != "") { ct.readTable(countfile, true, false);  }
		
		if (m->getControl_pressed()) { return 0; }
		
				
		if (listfile == "") {
			
			ofstream outSummary;
            map<string, string> variables; 
            variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(fastafile));
			string outputSummaryFile = getOutputFileName("summary", variables);
			util.openOutputFile(outputSummaryFile, outSummary);
			outSummary.setf(ios::fixed, ios::floatfield); outSummary.setf(ios::showpoint);
			outputNames.push_back(outputSummaryFile); outputTypes["summary"].push_back(outputSummaryFile);
			
			outSummary << "PositioninAlignment\tA\tT\tG\tC\tGap\tNumberofSeqs\tConsensusBase" << endl;
			
			ofstream outFasta;
			string outputFastaFile = getOutputFileName("fasta", variables);
			util.openOutputFile(outputFastaFile, outFasta);
			outputNames.push_back(outputFastaFile); outputTypes["fasta"].push_back(outputFastaFile);
        
			vector< vector<float> > percentages; percentages.resize(5);
			for (int j = 0; j < percentages.size(); j++) { percentages[j].resize(seqLength, 0.0); }
			
			string consSeq = "";
            int thisCount;
			//get counts
			for (int j = 0; j < seqLength; j++) {
				
				if (m->getControl_pressed()) { outSummary.close(); outFasta.close(); for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); } return 0; }
				
				vector<int> counts; counts.resize(5, 0); //A,T,G,C,Gap
				int numDots = 0;
				thisCount = 0;
				for (map<string, string>::iterator it = fastaMap.begin(); it != fastaMap.end(); it++) {
					
                    string thisSeq = it->second;
                    int size = 0;
                    
                    if (countfile != "") { size = ct.getNumSeqs(it->first); }
                    else {
                        map<string, int>::iterator itCount = nameFileMap.find(it->first);
                        if (itCount != nameFileMap.end()) {
                            size = itCount->second;
                        }else { m->mothurOut("[ERROR]: file mismatch, aborting.\n"); m->setControl_pressed(true); break; }
                    }
                    
                    for (int k = 0; k < size; k++) {
                        if (thisSeq[j] == '.') { numDots++; }
                        
                        char base = toupper(thisSeq[j]);
                        if (base == 'A') { counts[0]++; }
                        else if (base == 'T') { counts[1]++; }
                        else if (base == 'G') { counts[2]++; }
                        else if (base == 'C') { counts[3]++; }
                        else { counts[4]++; }
                        thisCount++;
                    }
				}
				
				char conBase = '.';
				if (numDots != thisCount) { conBase = getBase(counts, thisCount); }
				
				consSeq += conBase;
				
				percentages[0][j] = counts[0] / (float) thisCount;
				percentages[1][j] = counts[1] / (float) thisCount;
				percentages[2][j] = counts[2] / (float) thisCount;
				percentages[3][j] = counts[3] / (float) thisCount;
				percentages[4][j] = counts[4] / (float) thisCount;
			}
			
			for (int j = 0; j < seqLength; j++) { 
				outSummary << (j+1) << '\t' << percentages[0][j] << '\t'<< percentages[1][j] << '\t'<< percentages[2][j] << '\t' << percentages[3][j] << '\t' << percentages[4][j] << '\t' << thisCount << '\t' << consSeq[j] << endl;
			}
			
				
			outFasta << ">conseq" << endl << consSeq << endl;
			
			outSummary.close(); outFasta.close();
            
			if (m->getControl_pressed()) {  for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); } return 0; }
		
		}else {
			
            
			InputData* input = new InputData(listfile, "list", nullVector);
			ListVector* list = input->getListVector();
			
			string lastLabel = list->getLabel();
			set<string> processedLabels;
			set<string> userLabels = labels;

			//as long as you are not at the end of the file or done wih the lines you want
			while((list != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
				
				if (m->getControl_pressed()) {  for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); } delete list; delete input;  return 0;  }
				
				if(allLines == 1 || labels.count(list->getLabel()) == 1){			
					
					m->mothurOut(list->getLabel()); m->mothurOutEndLine();
					
					processList(list);
					
					processedLabels.insert(list->getLabel());
					userLabels.erase(list->getLabel());
				}
				
				if ((util.anyLabelsToProcess(list->getLabel(), userLabels, "") ) && (processedLabels.count(lastLabel) != 1)) {
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
			
			
			if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); } if (list != NULL) { delete list; } delete input; return 0;  }
			
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
			if (needToRun )  {
				if (list != NULL) { delete list; }
				
				list = input->getListVector(lastLabel);
				
				m->mothurOut(list->getLabel()); m->mothurOutEndLine();
				
				processList(list);
				
				delete list; list = NULL;
			}
			
			if (list != NULL) { delete list; }
			delete input;
		}
		
        m->mothurOut("It took " + toString(time(NULL) - start) + " secs to find the consensus sequences.");
        
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
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
        map<string, string> variables; 
        variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(fastafile));
        variables["[tag]"] = list->getLabel();
		string outputSummaryFile = getOutputFileName("summary", variables);
		util.openOutputFile(outputSummaryFile, outSummary);
		outSummary.setf(ios::fixed, ios::floatfield); outSummary.setf(ios::showpoint);
		outputNames.push_back(outputSummaryFile); outputTypes["summary"].push_back(outputSummaryFile);
		
		ofstream outName;
		string outputNameFile = getOutputFileName("name",variables);
		util.openOutputFile(outputNameFile, outName);
		outputNames.push_back(outputNameFile); outputTypes["name"].push_back(outputNameFile);
		
		ofstream outFasta;
		string outputFastaFile = getOutputFileName("fasta",variables);
		util.openOutputFile(outputFastaFile, outFasta);
		outputNames.push_back(outputFastaFile); outputTypes["fasta"].push_back(outputFastaFile);
		
		outSummary << "OTU#\tPositioninAlignment\tA\tT\tG\tC\tGap\tNumberofSeqs\tConsensusBase" << endl;
		
        string snumBins = toString(list->getNumBins());
        vector<string> binLabels = list->getLabels();
		for (int i = 0; i < list->getNumBins(); i++) {
			
			if (m->getControl_pressed()) { outSummary.close(); outName.close(); outFasta.close(); return 0; }
			
			string bin = list->get(i);
			string consSeq = getConsSeq(bin, outSummary, i);
			
			outFasta << ">" << binLabels[i] << endl << consSeq << endl;
			outName << binLabels[i] << '\t' << binLabels[i] << "," << bin << endl;
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
string ConsensusSeqsCommand::getConsSeq(string bin, ofstream& outSummary, int binNumber){
	try{
		
		string consSeq = "";
		bool error = false;
        int totalSize=0;
				
		vector<string> binNames;
		util.splitAtComma(bin, binNames);
        
        vector< vector<float> > percentages; percentages.resize(5);
		for (int j = 0; j < percentages.size(); j++) { percentages[j].resize(seqLength, 0.0); }

        if (countfile != "") {
            //get counts
            for (int j = 0; j < seqLength; j++) {
                
                if (m->getControl_pressed()) { return consSeq; }
                
                vector<int> counts; counts.resize(5, 0); //A,T,G,C,Gap
                int numDots = 0;
                totalSize = 0;
                 for (int i = 0; i < binNames.size(); i++) {
                     if (m->getControl_pressed()) { return consSeq; }
                     
                     string thisSeq = "";
                     map<string, string>::iterator itFasta = fastaMap.find(binNames[i]);
                     if (itFasta != fastaMap.end()) {
                         thisSeq = itFasta->second;
                     }else { m->mothurOut("[ERROR]: " + binNames[i] + " is not in your fasta file, please correct."); m->mothurOutEndLine(); m->setControl_pressed(true); }
                     
                     int size = ct.getNumSeqs(binNames[i]);
                     if (size != 0) {
                         for (int k = 0; k < size; k++) {
                             if (thisSeq[j] == '.') { numDots++; }
                             
                             char base = toupper(thisSeq[j]);
                             if (base == 'A') { counts[0]++; }
                             else if (base == 'T') { counts[1]++; }
                             else if (base == 'G') { counts[2]++; }
                             else if (base == 'C') { counts[3]++; }
                             else { counts[4]++; }
                             totalSize++;
                         }
                     }else { m->mothurOut("[ERROR]: " + binNames[i] + " is not in your count file, please correct."); m->mothurOutEndLine(); m->setControl_pressed(true); }
                 }
                char conBase = '.';
                if (numDots != totalSize) { conBase = getBase(counts, totalSize); }
                
                consSeq += conBase;
                
                percentages[0][j] = counts[0] / (float) totalSize;
                percentages[1][j] = counts[1] / (float) totalSize;
                percentages[2][j] = counts[2] / (float) totalSize;
                percentages[3][j] = counts[3] / (float) totalSize;
                percentages[4][j] = counts[4] / (float) totalSize;
            }

        }else {
		
            //get sequence strings for each name in the bin
            vector<string> seqs;
            for (int i = 0; i < binNames.size(); i++) {
                
                map<string, string>::iterator it;
                it = nameMap.find(binNames[i]);
                if (it == nameMap.end()) { 
                    if (namefile == "") { m->mothurOut("[ERROR]: " + binNames[i] + " is not in your fasta file, please correct."); m->mothurOutEndLine(); error = true; }
                    else { m->mothurOut("[ERROR]: " + binNames[i] + " is not in your fasta or name file, please correct."); m->mothurOutEndLine(); error = true; }
                    break;
                }else {
                    //add sequence string to seqs vector to process below
                    map<string, string>::iterator itFasta = fastaMap.find(it->second);
                    
                    if (itFasta != fastaMap.end()) {
                        string seq = itFasta->second;
                        seqs.push_back(seq);
                    }else { m->mothurOut("[ERROR]: file mismatch, aborting. \n"); }
                }
            }
            
            if (error) { m->setControl_pressed(true); return consSeq; }
            totalSize = seqs.size();
            //get counts
            for (int j = 0; j < seqLength; j++) {
                
                if (m->getControl_pressed()) { return consSeq; }
                
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
		}
        
        
        
		for (int j = 0; j < seqLength; j++) { 
			outSummary << (binNumber + 1) << '\t' << (j+1) << '\t' << percentages[0][j] << '\t'<< percentages[1][j] << '\t'<< percentages[2][j] << '\t' << percentages[3][j] << '\t' << percentages[4][j] << '\t' << totalSize << '\t' << consSeq[j] << endl;
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
        
		for (int i = 0; i < counts.size(); i++) {
            float countPercentage = counts[i] / (float) size;
			if (countPercentage < percentage) { counts[i] = 0; }
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
		util.openInputFile(fastafile, in);
		seqLength = 0;
        
		while (!in.eof()) {
			
			if (m->getControl_pressed()) { break; }
			
			Sequence seq(in); util.gobble(in);
			string name = seq.getName();
			
			if (name != "") {
				fastaMap[name] = seq.getAligned();
				nameMap[name] = name; //set nameMap incase no names file
				nameFileMap[name] = 1;
                
                if (seqLength == 0) { seqLength = seq.getAligned().length(); }
				else if (seqLength != seq.getAligned().length()) { m->mothurOut("[ERROR]: sequence are not the same length, please correct."); m->mothurOutEndLine(); m->setControl_pressed(true); break; }
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
         map<string, string> temp;
         map<string, string>::iterator it;
         bool error = false;
         
         util.readNames(namefile, temp); //use central buffered read
         
         for (map<string, string>::iterator itTemp = temp.begin(); itTemp != temp.end(); itTemp++) {
             string thisname, repnames;
             thisname = itTemp->first;
             repnames = itTemp->second;
             
             it = nameMap.find(thisname);
			 if (it != nameMap.end()) { //then this sequence was in the fastafile
				 nameFileMap[thisname] = util.getNumNames(repnames);	//for later when outputting the new namesFile if the list file is unique
                 
				 vector<string> splitRepNames;
				 util.splitAtComma(repnames, splitRepNames);
				 
				 for (int i = 0; i < splitRepNames.size(); i++) { nameMap[splitRepNames[i]] = thisname; }
				 
			 }else{	m->mothurOut("[ERROR]: " + thisname + " is not in the fasta file, please correct."); m->mothurOutEndLine(); error = true; }
         }
         
		 if (error) { m->setControl_pressed(true); }
 
		 return 0;
 
	}
	 catch(exception& e) {
		 m->errorOut(e, "ConsensusSeqsCommand", "readNames");
		 exit(1);
	 }
 }
 
//***************************************************************************************************************


