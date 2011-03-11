/*
 *  trimflowscommand.cpp
 *  Mothur
 *
 *  Created by Pat Schloss on 12/22/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "trimflowscommand.h"
#include "needlemanoverlap.hpp"

//**********************************************************************************************************************

vector<string> TrimFlowsCommand::getValidParameters(){	
	try {
		string Array[] =  {"flow", "maxflows", "minflows",
			"fasta", "minlength", "maxlength", "maxhomop", "signal", "noise"
			"oligos", "pdiffs", "bdiffs", "tdiffs",  "order",
			"allfiles", "processors",
			"outputdir","inputdir"
		
		};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "TrimFlowsCommand", "getValidParameters");
		exit(1);
	}
}

//**********************************************************************************************************************

vector<string> TrimFlowsCommand::getRequiredParameters(){	
	try {
		string Array[] =  {"flow"};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "TrimFlowsCommand", "getRequiredParameters");
		exit(1);
	}
}

//**********************************************************************************************************************

vector<string> TrimFlowsCommand::getRequiredFiles(){	
	try {
		vector<string> myArray;
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "TrimFlowsCommand", "getRequiredFiles");
		exit(1);
	}
}

//**********************************************************************************************************************

TrimFlowsCommand::TrimFlowsCommand(){	
	try {
		abort = true; calledHelp = true; 
		vector<string> tempOutNames;
		outputTypes["flow"] = tempOutNames;
		outputTypes["fasta"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "TrimFlowsCommand", "TrimFlowsCommand");
		exit(1);
	}
}

//***************************************************************************************************************

TrimFlowsCommand::~TrimFlowsCommand(){	/*	do nothing	*/	}

//***************************************************************************************************************

void TrimFlowsCommand::help(){
	try {
		m->mothurOut("The trim.flows command reads a flowgram file and creates .....\n");
		m->mothurOut("Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFasta).\n");
		m->mothurOut("For more details please check out the wiki http://www.mothur.org/wiki/Trim.flows.\n\n");
		
	}
	catch(exception& e) {
		m->errorOut(e, "TrimFlowsCommand", "help");
		exit(1);
	}
}

//**********************************************************************************************************************

TrimFlowsCommand::TrimFlowsCommand(string option)  {
	try {
		
		abort = false; calledHelp = false;   
		comboStarts = 0;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		
		else {
			//valid paramters for this command
			string AlignArray[] =  {"flow", "maxflows", "minflows",
				"fasta", "minlength", "maxlength", "maxhomop", "signal", "noise",
				"oligos", "pdiffs", "bdiffs", "tdiffs", "order",
				"allfiles", "processors",
		
				//			"group",
				"outputdir","inputdir"
				
			};
			
			vector<string> myArray (AlignArray, AlignArray+(sizeof(AlignArray)/sizeof(string)));
			
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
			outputTypes["flow"] = tempOutNames;
			outputTypes["fasta"] = tempOutNames;
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("flow");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["flow"] = inputDir + it->second;		}
				}
				
				it = parameters.find("oligos");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["oligos"] = inputDir + it->second;		}
				}
				
				
//				it = parameters.find("group");
//				//user has given a template file
//				if(it != parameters.end()){ 
//					path = m->hasPath(it->second);
//					//if the user has not given a path then, add inputdir. else leave path alone.
//					if (path == "") {	parameters["group"] = inputDir + it->second;		}
//				}
			}
			
			
			//check for required parameters
			flowFileName = validParameter.validFile(parameters, "flow", true);
			if (flowFileName == "not found") { m->mothurOut("flow is a required parameter for the trim.flows command."); m->mothurOutEndLine(); abort = true; }
			else if (flowFileName == "not open") { abort = true; }	
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += m->hasPath(flowFileName); //if user entered a file with a path then preserve it	
			}
			
			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			
			string temp;
			temp = validParameter.validFile(parameters, "minflows", false);	if (temp == "not found") { temp = "360"; }
			convert(temp, minFlows);  

			temp = validParameter.validFile(parameters, "maxflows", false);	if (temp == "not found") { temp = "720"; }
			convert(temp, maxFlows);  
			
			
			temp = validParameter.validFile(parameters, "oligos", true);
			if (temp == "not found")	{	oligoFileName = "";		}
			else if(temp == "not open")	{	abort = true;			} 
			else						{	oligoFileName = temp;	}
			
			temp = validParameter.validFile(parameters, "fasta", false);		if (temp == "not found"){	fasta = 0;		}
			else if(m->isTrue(temp))	{	fasta = 1;	}
			
			temp = validParameter.validFile(parameters, "maxhomop", false);		if (temp == "not found"){	temp = "9";		}
			convert(temp, maxHomoP);  

			temp = validParameter.validFile(parameters, "signal", false);		if (temp == "not found"){	temp = "0.50";	}
			convert(temp, signal);  

			temp = validParameter.validFile(parameters, "noise", false);		if (temp == "not found"){	temp = "0.70";	}
			convert(temp, noise);  

			temp = validParameter.validFile(parameters, "minlength", false);	if (temp == "not found"){	temp = "0";		}
			convert(temp, minLength); 
			
			temp = validParameter.validFile(parameters, "maxlength", false);	if (temp == "not found"){	temp = "0";		}
			convert(temp, maxLength);
			
			temp = validParameter.validFile(parameters, "bdiffs", false);		if (temp == "not found"){	temp = "0";		}
			convert(temp, bdiffs);
			
			temp = validParameter.validFile(parameters, "pdiffs", false);		if (temp == "not found"){	temp = "0";		}
			convert(temp, pdiffs);
			
			temp = validParameter.validFile(parameters, "tdiffs", false);
			if (temp == "not found"){ int tempTotal = pdiffs + bdiffs;  temp = toString(tempTotal); }
			convert(temp, tdiffs);
			if(tdiffs == 0){	tdiffs = bdiffs + pdiffs;	}
			
			temp = validParameter.validFile(parameters, "allfiles", false);		if (temp == "not found"){ temp = "T";		}
			allFiles = m->isTrue(temp);
			
			temp = validParameter.validFile(parameters, "processors", false);	if (temp == "not found"){ temp = "1";		}
			convert(temp, processors); 
	
			flowOrder = validParameter.validFile(parameters, "order", false);
			if (flowOrder == "not found"){ flowOrder = "TACG";		}
			else if(flowOrder.length() != 4){
				m->mothurOut("The value of the order option must be four bases long\n");
			}

			if(oligoFileName == ""){	allFiles = 0;		}

			numFPrimers = 0;
			numRPrimers = 0;
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "TrimFlowsCommand", "TrimSeqsCommand");
		exit(1);
	}
}

//***************************************************************************************************************

int TrimFlowsCommand::execute(){
	try{
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}

		string trimFlowFileName = outputDir + m->getRootName(m->getSimpleName(flowFileName)) + "trim.flow";
		outputNames.push_back(trimFlowFileName); outputTypes["flow"].push_back(trimFlowFileName);
		
		string scrapFlowFileName = outputDir + m->getRootName(m->getSimpleName(flowFileName)) + "scrap.flow";
		outputNames.push_back(scrapFlowFileName); outputTypes["flow"].push_back(scrapFlowFileName);

		string fastaFileName = outputDir + m->getRootName(m->getSimpleName(flowFileName)) + "flow.fasta";
		if(fasta){
			outputNames.push_back(fastaFileName); outputTypes["fasta"].push_back(fastaFileName);
		}
		
		vector<unsigned long int> flowFilePos = getFlowFileBreaks();
		for (int i = 0; i < (flowFilePos.size()-1); i++) {
			lines.push_back(new linePair(flowFilePos[i], flowFilePos[(i+1)]));
		}	

		vector<vector<string> > barcodePrimerComboFileNames;
		if(oligoFileName != ""){
			getOligos(barcodePrimerComboFileNames);	
		}
		
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
		if(processors == 1){
			driverCreateTrim(flowFileName, trimFlowFileName, scrapFlowFileName, fastaFileName, barcodePrimerComboFileNames, lines[0]);
		}else{
			createProcessesCreateTrim(flowFileName, trimFlowFileName, scrapFlowFileName, fastaFileName, barcodePrimerComboFileNames); 
		}	
#else
		driverCreateTrim(flowFileName, trimFlowFileName, scrapFlowFileName, fastaFileName, barcodePrimerComboFileNames, lines[0]);
#endif
		
		if (m->control_pressed) {  return 0; }			
		
		string flowFilesFileName;
		ofstream output;
		
		if(allFiles){
			
			flowFilesFileName = outputDir + m->getRootName(m->getSimpleName(flowFileName)) + "flow.files";
			m->openOutputFile(flowFilesFileName, output);

			for(int i=0;i<barcodePrimerComboFileNames.size();i++){
				for(int j=0;j<barcodePrimerComboFileNames[0].size();j++){
					
					FILE * pFile;
					unsigned long int size;
					
					//get num bytes in file
					pFile = fopen (barcodePrimerComboFileNames[i][j].c_str(),"rb");
					if (pFile==NULL) perror ("Error opening file");
					else{
						fseek (pFile, 0, SEEK_END);
						size=ftell (pFile);
						fclose (pFile);
					}

					if(size < 10){
						remove(barcodePrimerComboFileNames[i][j].c_str());
					}
					else{
						output << barcodePrimerComboFileNames[i][j] << endl;
					}
				}
			}
			output.close();
		}
		else{
			flowFilesFileName = outputDir + m->getRootName(m->getSimpleName(flowFileName)) + "flow.files";
			m->openOutputFile(flowFilesFileName, output);
			
			output << trimFlowFileName << endl;
			
			output.close();
		}
		outputTypes["flow.files"].push_back(flowFilesFileName);
		outputNames.push_back(flowFileName);
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
		
		return 0;	
	}
	catch(exception& e) {
		m->errorOut(e, "TrimSeqsCommand", "execute");
		exit(1);
	}
}

//***************************************************************************************************************

int TrimFlowsCommand::driverCreateTrim(string flowFileName, string trimFlowFileName, string scrapFlowFileName, string fastaFileName, vector<vector<string> > barcodePrimerComboFileNames, linePair* line){
	
	try {
		
		ofstream trimFlowFile;
		m->openOutputFile(trimFlowFileName, trimFlowFile);
		trimFlowFile.setf(ios::fixed, ios::floatfield); trimFlowFile.setf(ios::showpoint);

		ofstream scrapFlowFile;
		m->openOutputFile(scrapFlowFileName, scrapFlowFile);
		scrapFlowFile.setf(ios::fixed, ios::floatfield); scrapFlowFile.setf(ios::showpoint);

		if(line->start == 4){
			scrapFlowFile << maxFlows << endl;
			trimFlowFile << maxFlows << endl;
			if(allFiles){
				for(int i=0;i<barcodePrimerComboFileNames.size();i++){
					for(int j=0;j<barcodePrimerComboFileNames[0].size();j++){
						//				barcodePrimerComboFileNames[i][j] += toString(getpid()) + ".temp";
						ofstream temp;
						m->openOutputFile(barcodePrimerComboFileNames[i][j], temp);
						temp << maxFlows << endl;
						temp.close();
					}
				}			
			}
		}
		
		FlowData flowData(numFlows, signal, noise, maxHomoP, flowOrder);
		
		ofstream fastaFile;
		if(fasta){	m->openOutputFile(fastaFileName, fastaFile);	}
		
		ifstream flowFile;
		m->openInputFile(flowFileName, flowFile);
		
		flowFile.seekg(line->start);

		int count = 0;
		bool moreSeqs = 1;
			
		while(moreSeqs) {
			
			int success = 1;
			int currentSeqDiffs = 0;
			string trashCode = "";
			
			flowData.getNext(flowFile);
			flowData.capFlows(maxFlows);	
			
			Sequence currSeq = flowData.getSequence();
			if(!flowData.hasMinFlows(minFlows)){	//screen to see if sequence is of a minimum number of flows
				success = 0;
				trashCode += 'l';
			}
			
			if(minLength > 0 || maxLength > 0){	//screen to see if sequence is above and below a specific number of bases
				int seqLength = currSeq.getNumBases();
				if(seqLength < minLength || seqLength > maxLength){
					success = 0;
					trashCode += 'l';
				}
			}
			
			int primerIndex = 0;
			int barcodeIndex = 0;
			
			if(barcodes.size() != 0){
				success = stripBarcode(currSeq, barcodeIndex);
				if(success > bdiffs)		{	trashCode += 'b';	}
				else{ currentSeqDiffs += success;  }
			}
			
			if(numFPrimers != 0){
				success = stripForward(currSeq, primerIndex);
				if(success > pdiffs)		{	trashCode += 'f';	}
				else{ currentSeqDiffs += success;  }
			}
			
			if (currentSeqDiffs > tdiffs)	{	trashCode += 't';   }
			
			if(numRPrimers != 0){
				success = stripReverse(currSeq);
				if(!success)				{	trashCode += 'r';	}
			}

			if(trashCode.length() == 0){
							
				flowData.printFlows(trimFlowFile);
			
				if(fasta)	{	currSeq.printSequence(fastaFile);	}
				
				if(allFiles){
					ofstream output;
					m->openOutputFileAppend(barcodePrimerComboFileNames[barcodeIndex][primerIndex], output);
					output.setf(ios::fixed, ios::floatfield); trimFlowFile.setf(ios::showpoint);
					
					flowData.printFlows(output);
					output.close();
				}				
			}
			else{
				flowData.printFlows(scrapFlowFile, trashCode);
			}
				
			count++;
						
			//report progress
			if((count) % 10000 == 0){	m->mothurOut(toString(count)); m->mothurOutEndLine();		}

#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
			unsigned long int pos = flowFile.tellg();

			if ((pos == -1) || (pos >= line->end)) { break; }
#else
			if (flowFile.eof()) { break; }
#endif
			
		}
		//report progress
		if((count) % 10000 != 0){	m->mothurOut(toString(count)); m->mothurOutEndLine();		}
		
		trimFlowFile.close();
		scrapFlowFile.close();
		flowFile.close();
		if(fasta){	fastaFile.close();	}
		
		return count;
	}
	catch(exception& e) {
		m->errorOut(e, "TrimSeqsCommand", "driverCreateTrim");
		exit(1);
	}
}

//***************************************************************************************************************

void TrimFlowsCommand::getOligos(vector<vector<string> >& outFlowFileNames){
	try {
		ifstream oligosFile;
		m->openInputFile(oligoFileName, oligosFile);
		
		string type, oligo, group;

		int indexPrimer = 0;
		int indexBarcode = 0;
		
		while(!oligosFile.eof()){
		
			oligosFile >> type; m->gobble(oligosFile);	//get the first column value of the row - is it a comment or a feature we are interested in?

			if(type[0] == '#'){	//igore the line because there's a comment
				while (!oligosFile.eof())	{	char c = oligosFile.get(); if (c == 10 || c == 13){	break;	}	} // get rest of line if there's any crap there
			}
			else{				//there's a feature we're interested in

				for(int i=0;i<type.length();i++){	type[i] = toupper(type[i]);  }					//make type case insensitive

				oligosFile >> oligo;	//get the DNA sequence for the feature

				for(int i=0;i<oligo.length();i++){	//make type case insensitive and change any U's to T's
					oligo[i] = toupper(oligo[i]);
					if(oligo[i] == 'U')	{	oligo[i] = 'T';	}
				}

				if(type == "FORWARD"){	//if the feature is a forward primer...
					group = "";

					while (!oligosFile.eof())	{	// get rest of line in case there is a primer name = will have the name of the primer
						char c = oligosFile.get(); 
						if (c == 10 || c == 13){	break;	}
						else if (c == 32 || c == 9){;} //space or tab
						else { 	group += c;  }
					} 

					//have we seen this primer already?
					map<string, int>::iterator itPrimer = primers.find(oligo);
					if (itPrimer != primers.end()) { m->mothurOut("primer " + oligo + " is in your oligos file already."); m->mothurOutEndLine();  }

					primers[oligo]=indexPrimer; indexPrimer++;
					primerNameVector.push_back(group);

				}
				else if(type == "REVERSE"){
					Sequence oligoRC("reverse", oligo);
					oligoRC.reverseComplement();
					revPrimer.push_back(oligoRC.getUnaligned());
				}
				else if(type == "BARCODE"){
					oligosFile >> group;

					//check for repeat barcodes
					map<string, int>::iterator itBar = barcodes.find(oligo);
					if (itBar != barcodes.end()) { m->mothurOut("barcode " + oligo + " is in your oligos file already."); m->mothurOutEndLine();  }

					barcodes[oligo]=indexBarcode; indexBarcode++;
					barcodeNameVector.push_back(group);
				}
				else{
					m->mothurOut(type + " is not recognized as a valid type. Choices are forward, reverse, and barcode. Ignoring " + oligo + "."); m->mothurOutEndLine();  
				}
			}

			m->gobble(oligosFile);
		}
		oligosFile.close();
		
		if(barcodeNameVector.size() == 0 && primerNameVector[0] == ""){	allFiles = 0;	}
		
		//add in potential combos
		if(barcodeNameVector.size() == 0){
			barcodes[""] = 0;
			barcodeNameVector.push_back("");			
		}
		
		if(primerNameVector.size() == 0){
			primers[""] = 0;
			primerNameVector.push_back("");			
		}
		
		
		outFlowFileNames.resize(barcodeNameVector.size());
		for(int i=0;i<outFlowFileNames.size();i++){
			outFlowFileNames[i].assign(primerNameVector.size(), "");
		}
		
		if(allFiles){

			for(map<string, int>::iterator itBar = barcodes.begin();itBar != barcodes.end();itBar++){
				for(map<string, int>::iterator itPrimer = primers.begin();itPrimer != primers.end(); itPrimer++){

					string primerName = primerNameVector[itPrimer->second];
					string barcodeName = barcodeNameVector[itBar->second];
										
					string comboGroupName = "";
					string fileName = "";
					
					if(primerName == ""){
						comboGroupName = barcodeNameVector[itBar->second];
						fileName = outputDir + m->getRootName(m->getSimpleName(flowFileName)) + comboGroupName + ".flow";
					}
					else{
						if(barcodeName == ""){
							comboGroupName = primerNameVector[itPrimer->second];
						}
						else{
							comboGroupName = barcodeNameVector[itBar->second] + "." + primerNameVector[itPrimer->second];
						}
						fileName = outputDir + m->getRootName(m->getSimpleName(flowFileName)) + comboGroupName + ".flow";
					}
					
					outputNames.push_back(fileName);
					outputTypes["flow"].push_back(fileName);
					outFlowFileNames[itBar->second][itPrimer->second] = fileName;
					
					ofstream temp;
					m->openOutputFile(fileName, temp);
					temp.close();
				}
			}
		}
		
		numFPrimers = primers.size();
		numRPrimers = revPrimer.size();
		
	}
	catch(exception& e) {
		m->errorOut(e, "TrimSeqsCommand", "getOligos");
		exit(1);
	}
}

//***************************************************************************************************************

int TrimFlowsCommand::stripBarcode(Sequence& seq, int& group){
	try {
		
		string rawSequence = seq.getUnaligned();
		int success = bdiffs + 1;	//guilty until proven innocent
		
		//can you find the barcode
		for(map<string,int>::iterator it=barcodes.begin();it!=barcodes.end();it++){
			string oligo = it->first;
			if(rawSequence.length() < oligo.length()){	//let's just assume that the barcodes are the same length
				success = bdiffs + 10;					//if the sequence is shorter than the barcode then bail out
				break;	
			}
			
			if(compareDNASeq(oligo, rawSequence.substr(0,oligo.length()))){
				group = it->second;
				seq.setUnaligned(rawSequence.substr(oligo.length()));
				
				success = 0;
				break;
			}
		}
		
		//if you found the barcode or if you don't want to allow for diffs
		if ((bdiffs == 0) || (success == 0)) { return success;  }
		
		else { //try aligning and see if you can find it
			
			int maxLength = 0;
			
			Alignment* alignment;
			if (barcodes.size() > 0) {
				map<string,int>::iterator it=barcodes.begin();
				
				for(it;it!=barcodes.end();it++){
					if(it->first.length() > maxLength){
						maxLength = it->first.length();
					}
				}
				alignment = new NeedlemanOverlap(-1.0, 1.0, -1.0, (maxLength+bdiffs+1));  
				
			}else{ alignment = NULL; } 
			
			//can you find the barcode
			int minDiff = 1e6;
			int minCount = 1;
			int minGroup = -1;
			int minPos = 0;
			
			for(map<string,int>::iterator it=barcodes.begin();it!=barcodes.end();it++){
				string oligo = it->first;
				//				int length = oligo.length();
				
				if(rawSequence.length() < maxLength){	//let's just assume that the barcodes are the same length
					success = bdiffs + 10;
					break;
				}
				
				//use needleman to align first barcode.length()+numdiffs of sequence to each barcode
				alignment->align(oligo, rawSequence.substr(0,oligo.length()+bdiffs));
				oligo = alignment->getSeqAAln();
				string temp = alignment->getSeqBAln();
				
				int alnLength = oligo.length();
				
				for(int i=oligo.length()-1;i>=0;i--){
					if(oligo[i] != '-'){	alnLength = i+1;	break;	}
				}
				oligo = oligo.substr(0,alnLength);
				temp = temp.substr(0,alnLength);
				
				int numDiff = countDiffs(oligo, temp);
				
				if(numDiff < minDiff){
					minDiff = numDiff;
					minCount = 1;
					minGroup = it->second;
					minPos = 0;
					for(int i=0;i<alnLength;i++){
						if(temp[i] != '-'){
							minPos++;
						}
					}
				}
				else if(numDiff == minDiff){
					minCount++;
				}
				
			}
			
			if(minDiff > bdiffs)	{	success = minDiff;		}	//no good matches
			else if(minCount > 1)	{	success = bdiffs + 100;	}	//can't tell the difference between multiple barcodes
			else{													//use the best match
				group = minGroup;
				seq.setUnaligned(rawSequence.substr(minPos));
				success = minDiff;
			}
			
			if (alignment != NULL) {  delete alignment;  }
			
		}
		
		return success;
		
	}
	catch(exception& e) {
		m->errorOut(e, "TrimFlowsCommand", "stripBarcode");
		exit(1);
	}
	
}

//***************************************************************************************************************

int TrimFlowsCommand::stripForward(Sequence& seq, int& group){
	try {
		
		string rawSequence = seq.getUnaligned();
		int success = pdiffs + 1;	//guilty until proven innocent
		
		//can you find the primer
		for(map<string,int>::iterator it=primers.begin();it!=primers.end();it++){
			string oligo = it->first;
			if(rawSequence.length() < oligo.length()){	//let's just assume that the primers are the same length
				success = pdiffs + 10;					//if the sequence is shorter than the barcode then bail out
				break;	
			}
			
			if(compareDNASeq(oligo, rawSequence.substr(0,oligo.length()))){
				group = it->second;
				seq.setUnaligned(rawSequence.substr(oligo.length()));
				success = 0;
				break;
			}
		}
		
		//if you found the barcode or if you don't want to allow for diffs
		if ((pdiffs == 0) || (success == 0)) {	return success;  }
		
		else { //try aligning and see if you can find it
			
			int maxLength = 0;
			
			Alignment* alignment;
			if (primers.size() > 0) {
				map<string,int>::iterator it=primers.begin();
				
				for(it;it!=primers.end();it++){
					if(it->first.length() > maxLength){
						maxLength = it->first.length();
					}
				}
				alignment = new NeedlemanOverlap(-1.0, 1.0, -1.0, (maxLength+pdiffs+1));  
				
			}else{ alignment = NULL; } 
			
			//can you find the barcode
			int minDiff = 1e6;
			int minCount = 1;
			int minGroup = -1;
			int minPos = 0;
			
			for(map<string,int>::iterator it=primers.begin();it!=primers.end();it++){
				string oligo = it->first;
				//				int length = oligo.length();
				
				if(rawSequence.length() < maxLength){	
					success = pdiffs + 100;
					break;
				}
				
				//use needleman to align first barcode.length()+numdiffs of sequence to each barcode
				alignment->align(oligo, rawSequence.substr(0,oligo.length()+pdiffs));
				oligo = alignment->getSeqAAln();
				string temp = alignment->getSeqBAln();
				
				int alnLength = oligo.length();
				
				for(int i=oligo.length()-1;i>=0;i--){
					if(oligo[i] != '-'){	alnLength = i+1;	break;	}
				}
				oligo = oligo.substr(0,alnLength);
				temp = temp.substr(0,alnLength);
				
				int numDiff = countDiffs(oligo, temp);
				
				if(numDiff < minDiff){
					minDiff = numDiff;
					minCount = 1;
					minGroup = it->second;
					minPos = 0;
					for(int i=0;i<alnLength;i++){
						if(temp[i] != '-'){
							minPos++;
						}
					}
				}
				else if(numDiff == minDiff){
					minCount++;
				}
				
			}
			
			if(minDiff > pdiffs)	{	success = minDiff;		}	//no good matches
			else if(minCount > 1)	{	success = pdiffs + 10;	}	//can't tell the difference between multiple primers
			else{													//use the best match
				group = minGroup;
				seq.setUnaligned(rawSequence.substr(minPos));
				success = minDiff;
			}
			
			if (alignment != NULL) {  delete alignment;  }
			
		}
		
		return success;
		
	}
	catch(exception& e) {
		m->errorOut(e, "TrimFlowsCommand", "stripForward");
		exit(1);
	}
}

//***************************************************************************************************************

bool TrimFlowsCommand::stripReverse(Sequence& seq){
	try {

		string rawSequence = seq.getUnaligned();
		bool success = 0;	//guilty until proven innocent
		
		for(int i=0;i<numRPrimers;i++){
			string oligo = revPrimer[i];
			
			if(rawSequence.length() < oligo.length()){
				success = 0;
				break;
			}
			
			if(compareDNASeq(oligo, rawSequence.substr(rawSequence.length()-oligo.length(),oligo.length()))){
				seq.setUnaligned(rawSequence.substr(0,rawSequence.length()-oligo.length()));
				success = 1;
				break;
			}
		}	

		return success;
		
	}
	catch(exception& e) {
		m->errorOut(e, "TrimFlowsCommand", "stripReverse");
		exit(1);
	}
}


//***************************************************************************************************************

bool TrimFlowsCommand::compareDNASeq(string oligo, string seq){
	try {
		bool success = 1;
		int length = oligo.length();
		
		for(int i=0;i<length;i++){
			
			if(oligo[i] != seq[i]){
				if(oligo[i] == 'A' || oligo[i] == 'T' || oligo[i] == 'G' || oligo[i] == 'C')	{	success = 0; 	}
				else if((oligo[i] == 'N' || oligo[i] == 'I') && (seq[i] == 'N'))				{	success = 0;	}
				else if(oligo[i] == 'R' && (seq[i] != 'A' && seq[i] != 'G'))					{	success = 0;	}
				else if(oligo[i] == 'Y' && (seq[i] != 'C' && seq[i] != 'T'))					{	success = 0;	}
				else if(oligo[i] == 'M' && (seq[i] != 'C' && seq[i] != 'A'))					{	success = 0;	}
				else if(oligo[i] == 'K' && (seq[i] != 'T' && seq[i] != 'G'))					{	success = 0;	}
				else if(oligo[i] == 'W' && (seq[i] != 'T' && seq[i] != 'A'))					{	success = 0;	}
				else if(oligo[i] == 'S' && (seq[i] != 'C' && seq[i] != 'G'))					{	success = 0;	}
				else if(oligo[i] == 'B' && (seq[i] != 'C' && seq[i] != 'T' && seq[i] != 'G'))	{	success = 0;	}
				else if(oligo[i] == 'D' && (seq[i] != 'A' && seq[i] != 'T' && seq[i] != 'G'))	{	success = 0;	}
				else if(oligo[i] == 'H' && (seq[i] != 'A' && seq[i] != 'T' && seq[i] != 'C'))	{	success = 0;	}
				else if(oligo[i] == 'V' && (seq[i] != 'A' && seq[i] != 'C' && seq[i] != 'G'))	{	success = 0;	}			
				
				if(success == 0)	{	break;	 }
			}
			else{
				success = 1;
			}
		}
		
		return success;
	}
	catch(exception& e) {
		m->errorOut(e, "TrimFlowsCommand", "compareDNASeq");
		exit(1);
	}
	
}

//***************************************************************************************************************

int TrimFlowsCommand::countDiffs(string oligo, string seq){
	try {
		
		int length = oligo.length();
		int countDiffs = 0;
		
		for(int i=0;i<length;i++){
			
			if(oligo[i] != seq[i]){
				if(oligo[i] == 'A' || oligo[i] == 'T' || oligo[i] == 'G' || oligo[i] == 'C' || oligo[i] == '-' || oligo[i] == '.')	{	countDiffs++; 	}
				else if((oligo[i] == 'N' || oligo[i] == 'I') && (seq[i] == 'N'))				{	countDiffs++;	}
				else if(oligo[i] == 'R' && (seq[i] != 'A' && seq[i] != 'G'))					{	countDiffs++;	}
				else if(oligo[i] == 'Y' && (seq[i] != 'C' && seq[i] != 'T'))					{	countDiffs++;	}
				else if(oligo[i] == 'M' && (seq[i] != 'C' && seq[i] != 'A'))					{	countDiffs++;	}
				else if(oligo[i] == 'K' && (seq[i] != 'T' && seq[i] != 'G'))					{	countDiffs++;	}
				else if(oligo[i] == 'W' && (seq[i] != 'T' && seq[i] != 'A'))					{	countDiffs++;	}
				else if(oligo[i] == 'S' && (seq[i] != 'C' && seq[i] != 'G'))					{	countDiffs++;	}
				else if(oligo[i] == 'B' && (seq[i] != 'C' && seq[i] != 'T' && seq[i] != 'G'))	{	countDiffs++;	}
				else if(oligo[i] == 'D' && (seq[i] != 'A' && seq[i] != 'T' && seq[i] != 'G'))	{	countDiffs++;	}
				else if(oligo[i] == 'H' && (seq[i] != 'A' && seq[i] != 'T' && seq[i] != 'C'))	{	countDiffs++;	}
				else if(oligo[i] == 'V' && (seq[i] != 'A' && seq[i] != 'C' && seq[i] != 'G'))	{	countDiffs++;	}	
			}
			
		}
		
		return countDiffs;
	}
	catch(exception& e) {
		m->errorOut(e, "TrimFlowsCommand", "countDiffs");
		exit(1);
	}
	
}

/**************************************************************************************************/

vector<unsigned long int> TrimFlowsCommand::getFlowFileBreaks() {

	try{
			
		vector<unsigned long int> filePos;
		filePos.push_back(0);
					
		FILE * pFile;
		unsigned long int size;
		
		//get num bytes in file
		pFile = fopen (flowFileName.c_str(),"rb");
		if (pFile==NULL) perror ("Error opening file");
		else{
			fseek (pFile, 0, SEEK_END);
			size=ftell (pFile);
			fclose (pFile);
		}
				
		//estimate file breaks
		unsigned long int chunkSize = 0;
		chunkSize = size / processors;

		//file too small to divide by processors
		if (chunkSize == 0)  {  processors = 1;	filePos.push_back(size); return filePos;	}
		
		//for each process seekg to closest file break and search for next '>' char. make that the filebreak
		for (int i = 0; i < processors; i++) {
			unsigned long int spot = (i+1) * chunkSize;
			
			ifstream in;
			m->openInputFile(flowFileName, in);
			in.seekg(spot);
			
			string dummy = m->getline(in);
			
			//there was not another sequence before the end of the file
			unsigned long int sanityPos = in.tellg();
			
//			if (sanityPos == -1) {	break;  }
//			else {  filePos.push_back(newSpot);  }
			if (sanityPos == -1) {	break;  }
			else {  filePos.push_back(sanityPos);  }
			
			in.close();
		}
		
		//save end pos
		filePos.push_back(size);
		
		//sanity check filePos
		for (int i = 0; i < (filePos.size()-1); i++) {
			if (filePos[(i+1)] <= filePos[i]) {  filePos.erase(filePos.begin()+(i+1)); i--; }
		}

		ifstream in;
		m->openInputFile(flowFileName, in);
		in >> numFlows;
		m->gobble(in);
		unsigned long int spot = in.tellg();
		filePos[0] = spot;
		in.close();
		
		processors = (filePos.size() - 1);
		
		return filePos;	
	}
	catch(exception& e) {
		m->errorOut(e, "TrimSeqsCommand", "getFlowFileBreaks");
		exit(1);
	}
}

/**************************************************************************************************/

int TrimFlowsCommand::createProcessesCreateTrim(string flowFileName, string trimFlowFileName, string scrapFlowFileName, string fastaFileName, vector<vector<string> > barcodePrimerComboFileNames){

	try {
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
		int process = 1;
		int exitCommand = 1;
		processIDS.clear();
		
		//loop through and create all the processes you want
		while (process != processors) {
			int pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
				process++;
			}else if (pid == 0){
				
				vector<vector<string> > tempBarcodePrimerComboFileNames = barcodePrimerComboFileNames;
				if(allFiles){
					for(int i=0;i<tempBarcodePrimerComboFileNames.size();i++){
						for(int j=0;j<tempBarcodePrimerComboFileNames[0].size();j++){
							tempBarcodePrimerComboFileNames[i][j] += toString(getpid()) + ".temp";
							ofstream temp;
							m->openOutputFile(tempBarcodePrimerComboFileNames[i][j], temp);
							temp.close();
							
						}
					}
				}
				driverCreateTrim(flowFileName,
								 (trimFlowFileName + toString(getpid()) + ".temp"),
								 (scrapFlowFileName + toString(getpid()) + ".temp"),
								 (fastaFileName + toString(getpid()) + ".temp"),
								 tempBarcodePrimerComboFileNames, lines[process]);

				exit(0);
			}else { 
				m->mothurOut("[ERROR]: unable to spawn the necessary processes."); m->mothurOutEndLine(); 
				for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
				exit(0);
			}
		}
		
		//parent do my part
		ofstream temp;
		m->openOutputFile(trimFlowFileName, temp);
		temp.close();

		m->openOutputFile(scrapFlowFileName, temp);
		temp.close();
		
		if(fasta){
			m->openOutputFile(fastaFileName, temp);
			temp.close();
		}
		
		driverCreateTrim(flowFileName, trimFlowFileName, scrapFlowFileName, fastaFileName, barcodePrimerComboFileNames, lines[0]);

		//force parent to wait until all the processes are done
		for (int i=0;i<processIDS.size();i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
		
		//append files
		m->mothurOutEndLine();
		for(int i=0;i<processIDS.size();i++){
			
			m->mothurOut("Appending files from process " + toString(processIDS[i])); m->mothurOutEndLine();
			
			m->appendFiles((trimFlowFileName + toString(processIDS[i]) + ".temp"), trimFlowFileName);
			remove((trimFlowFileName + toString(processIDS[i]) + ".temp").c_str());
//			m->mothurOut("\tDone with trim.flow file"); m->mothurOutEndLine();

			m->appendFiles((scrapFlowFileName + toString(processIDS[i]) + ".temp"), scrapFlowFileName);
			remove((scrapFlowFileName + toString(processIDS[i]) + ".temp").c_str());
//			m->mothurOut("\tDone with scrap.flow file"); m->mothurOutEndLine();

			if(fasta){
				m->appendFiles((fastaFileName + toString(processIDS[i]) + ".temp"), fastaFileName);
				remove((fastaFileName + toString(processIDS[i]) + ".temp").c_str());
//				m->mothurOut("\tDone with flow.fasta file"); m->mothurOutEndLine();
			}
			if(allFiles){						
				for (int j = 0; j < barcodePrimerComboFileNames.size(); j++) {
					for (int k = 0; k < barcodePrimerComboFileNames[0].size(); k++) {
						m->appendFiles((barcodePrimerComboFileNames[j][k] + toString(processIDS[i]) + ".temp"), barcodePrimerComboFileNames[j][k]);
						remove((barcodePrimerComboFileNames[j][k] + toString(processIDS[i]) + ".temp").c_str());
					}
				}
			}
		}
		
		return exitCommand;
#endif		
	}
	catch(exception& e) {
		m->errorOut(e, "TrimFlowsCommand", "createProcessesCreateTrim");
		exit(1);
	}
}

//***************************************************************************************************************
