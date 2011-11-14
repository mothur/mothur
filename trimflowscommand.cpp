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
vector<string> TrimFlowsCommand::setParameters(){	
	try {
		CommandParameter pflow("flow", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(pflow);
		CommandParameter poligos("oligos", "InputTypes", "", "", "none", "none", "none",false,false); parameters.push_back(poligos);
		CommandParameter pmaxhomop("maxhomop", "Number", "", "9", "", "", "",false,false); parameters.push_back(pmaxhomop);
		CommandParameter pmaxflows("maxflows", "Number", "", "450", "", "", "",false,false); parameters.push_back(pmaxflows);
		CommandParameter pminflows("minflows", "Number", "", "450", "", "", "",false,false); parameters.push_back(pminflows);
		CommandParameter ppdiffs("pdiffs", "Number", "", "0", "", "", "",false,false); parameters.push_back(ppdiffs);
		CommandParameter pbdiffs("bdiffs", "Number", "", "0", "", "", "",false,false); parameters.push_back(pbdiffs);
		CommandParameter ptdiffs("tdiffs", "Number", "", "0", "", "", "",false,false); parameters.push_back(ptdiffs);
		CommandParameter pprocessors("processors", "Number", "", "1", "", "", "",false,false); parameters.push_back(pprocessors);
		CommandParameter psignal("signal", "Number", "", "0.50", "", "", "",false,false); parameters.push_back(psignal);
		CommandParameter pnoise("noise", "Number", "", "0.70", "", "", "",false,false); parameters.push_back(pnoise);
		CommandParameter pallfiles("allfiles", "Boolean", "", "t", "", "", "",false,false); parameters.push_back(pallfiles);
		CommandParameter porder("order", "String", "", "", "", "", "",false,false); parameters.push_back(porder);
		CommandParameter pfasta("fasta", "Boolean", "", "F", "", "", "",false,false); parameters.push_back(pfasta);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "TrimFlowsCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string TrimFlowsCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The trim.flows command reads a flowgram file and creates .....\n";
		helpString += "Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFasta).\n";
		helpString += "For more details please check out the wiki http://www.mothur.org/wiki/Trim.flows.\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "TrimFlowsCommand", "getHelpString");
		exit(1);
	}
}

//**********************************************************************************************************************

TrimFlowsCommand::TrimFlowsCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["flow"] = tempOutNames;
		outputTypes["fasta"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "TrimFlowsCommand", "TrimFlowsCommand");
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
				
			}
			
			
			//check for required parameters
			flowFileName = validParameter.validFile(parameters, "flow", true);
			if (flowFileName == "not found") { 
				flowFileName = m->getFlowFile(); 
				if (flowFileName != "") {  m->mothurOut("Using " + flowFileName + " as input file for the flow parameter."); m->mothurOutEndLine(); }
				else { 
					m->mothurOut("No valid current flow file. You must provide a flow file."); m->mothurOutEndLine(); 
					abort = true;
				} 
			}else if (flowFileName == "not open") { flowFileName = ""; abort = true; }	
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += m->hasPath(flowFileName); //if user entered a file with a path then preserve it	
			}
			
			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			
			string temp;
			temp = validParameter.validFile(parameters, "minflows", false);	if (temp == "not found") { temp = "450"; }
			convert(temp, minFlows);  

			temp = validParameter.validFile(parameters, "maxflows", false);	if (temp == "not found") { temp = "450"; }
			convert(temp, maxFlows);  
			
			
			temp = validParameter.validFile(parameters, "oligos", true);
			if (temp == "not found")	{	oligoFileName = "";		}
			else if(temp == "not open")	{	abort = true;			} 
			else						{	oligoFileName = temp;	m->setOligosFile(oligoFileName); }
			
			temp = validParameter.validFile(parameters, "fasta", false);		if (temp == "not found"){	fasta = 0;		}
			else if(m->isTrue(temp))	{	fasta = 1;	}
			
			temp = validParameter.validFile(parameters, "maxhomop", false);		if (temp == "not found"){	temp = "9";		}
			convert(temp, maxHomoP);  

			temp = validParameter.validFile(parameters, "signal", false);		if (temp == "not found"){	temp = "0.50";	}
			convert(temp, signal);  

			temp = validParameter.validFile(parameters, "noise", false);		if (temp == "not found"){	temp = "0.70";	}
			convert(temp, noise);  
	
			temp = validParameter.validFile(parameters, "bdiffs", false);		if (temp == "not found"){	temp = "0";		}
			convert(temp, bdiffs);
			
			temp = validParameter.validFile(parameters, "pdiffs", false);		if (temp == "not found"){	temp = "0";		}
			convert(temp, pdiffs);
			
			temp = validParameter.validFile(parameters, "tdiffs", false);
			if (temp == "not found"){ int tempTotal = pdiffs + bdiffs;  temp = toString(tempTotal); }
			convert(temp, tdiffs);
			if(tdiffs == 0){	tdiffs = bdiffs + pdiffs;	}
			
			temp = validParameter.validFile(parameters, "processors", false);	if (temp == "not found"){	temp = m->getProcessors();	}
			m->setProcessors(temp);
			convert(temp, processors);
	
			flowOrder = validParameter.validFile(parameters, "order", false);
			if (flowOrder == "not found"){ flowOrder = "TACG";		}
			else if(flowOrder.length() != 4){
				m->mothurOut("The value of the order option must be four bases long\n");
			}

			if(oligoFileName == "")	{	allFiles = 0;		}
			else					{	allFiles = 1;		}

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
		
		vector<unsigned long long> flowFilePos;
	#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
		flowFilePos = getFlowFileBreaks();
		for (int i = 0; i < (flowFilePos.size()-1); i++) {
			lines.push_back(new linePair(flowFilePos[i], flowFilePos[(i+1)]));
		}	
	#else
		ifstream in; m->openInputFile(flowFileName, in); in >> numFlows; in.close();
	///////////////////////////////////////// until I fix multiple processors for windows //////////////////	
		processors = 1;
	///////////////////////////////////////// until I fix multiple processors for windows //////////////////		
		if (processors == 1) {
			lines.push_back(new linePair(0, 1000));
		}else {
			int numFlowLines;
			flowFilePos = m->setFilePosEachLine(flowFileName, numFlowLines);
			flowFilePos.erase(flowFilePos.begin() + 1); numFlowLines--;
			
			//figure out how many sequences you have to process
			int numSeqsPerProcessor = numFlowLines / processors;
			cout << numSeqsPerProcessor << '\t' << numFlowLines << endl;
			for (int i = 0; i < processors; i++) {
				int startIndex =  i * numSeqsPerProcessor;
				if(i == (processors - 1)){	numSeqsPerProcessor = numFlowLines - i * numSeqsPerProcessor; 	}
				lines.push_back(new linePair(flowFilePos[startIndex], numSeqsPerProcessor));
				cout << flowFilePos[startIndex] << '\t' << numSeqsPerProcessor << endl;
			}
		}
	#endif
		
		vector<vector<string> > barcodePrimerComboFileNames;
		if(oligoFileName != ""){
			getOligos(barcodePrimerComboFileNames);	
		}
		
		if(processors == 1){
			driverCreateTrim(flowFileName, trimFlowFileName, scrapFlowFileName, fastaFileName, barcodePrimerComboFileNames, lines[0]);
		}else{
			createProcessesCreateTrim(flowFileName, trimFlowFileName, scrapFlowFileName, fastaFileName, barcodePrimerComboFileNames); 
		}	
		
		if (m->control_pressed) {  return 0; }			
		
		string flowFilesFileName;
		ofstream output;
		
		if(allFiles){
			set<string> namesAlreadyProcessed;
			flowFilesFileName = outputDir + m->getRootName(m->getSimpleName(flowFileName)) + "flow.files";
			m->openOutputFile(flowFilesFileName, output);

			for(int i=0;i<barcodePrimerComboFileNames.size();i++){
				for(int j=0;j<barcodePrimerComboFileNames[0].size();j++){
					if (namesAlreadyProcessed.count(barcodePrimerComboFileNames[i][j]) == 0) {
						FILE * pFile;
						unsigned long long size;
						
						//get num bytes in file
						pFile = fopen (barcodePrimerComboFileNames[i][j].c_str(),"rb");
						if (pFile==NULL) perror ("Error opening file");
						else{
							fseek (pFile, 0, SEEK_END);
							size=ftell(pFile);
							fclose (pFile);
						}
						
						if(size < 10){
							m->mothurRemove(barcodePrimerComboFileNames[i][j]);
						}
						else{
							output << barcodePrimerComboFileNames[i][j] << endl;
							outputNames.push_back(barcodePrimerComboFileNames[i][j]);
							outputTypes["flow"].push_back(barcodePrimerComboFileNames[i][j]);
						}
						namesAlreadyProcessed.insert(barcodePrimerComboFileNames[i][j]);
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
		outputNames.push_back(flowFilesFileName);
		
//		set fasta file as new current fastafile
//		string current = "";
//		itTypes = outputTypes.find("fasta");
//		if (itTypes != outputTypes.end()) {
//			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setFastaFile(current); }
//		}
		
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

int TrimFlowsCommand::driverCreateTrim(string flowFileName, string trimFlowFileName, string scrapFlowFileName, string fastaFileName, vector<vector<string> > thisBarcodePrimerComboFileNames, linePair* line){
	
	try {
		ofstream trimFlowFile;
		m->openOutputFile(trimFlowFileName, trimFlowFile);
		trimFlowFile.setf(ios::fixed, ios::floatfield); trimFlowFile.setf(ios::showpoint);

		ofstream scrapFlowFile;
		m->openOutputFile(scrapFlowFileName, scrapFlowFile);
		scrapFlowFile.setf(ios::fixed, ios::floatfield); scrapFlowFile.setf(ios::showpoint);
		
		ofstream fastaFile;
		if(fasta){	m->openOutputFile(fastaFileName, fastaFile);	}
		
		ifstream flowFile;
		m->openInputFile(flowFileName, flowFile);
		
		flowFile.seekg(line->start);
		
		if(line->start == 0){
			flowFile >> numFlows; m->gobble(flowFile);
			scrapFlowFile << maxFlows << endl;
			trimFlowFile << maxFlows << endl;
			if(allFiles){
				for(int i=0;i<thisBarcodePrimerComboFileNames.size();i++){
					for(int j=0;j<thisBarcodePrimerComboFileNames[0].size();j++){
						ofstream temp;
						m->openOutputFile(thisBarcodePrimerComboFileNames[i][j], temp);
						temp << maxFlows << endl;
						temp.close();
					}
				}			
			}
		}
		
		FlowData flowData(numFlows, signal, noise, maxHomoP, flowOrder);
		//cout << " driver flowdata address " <<  &flowData  << &flowFile << endl;	
		int count = 0;
		bool moreSeqs = 1;
		
		TrimOligos trimOligos(pdiffs, bdiffs, primers, barcodes, revPrimer);
		
		while(moreSeqs) {
			//cout << "driver " << count << endl;
			
	
			if (m->control_pressed) { break; }
			
			int success = 1;
			int currentSeqDiffs = 0;
			string trashCode = "";
			
			flowData.getNext(flowFile); 
			//cout << "driver good bit " << flowFile.good() << endl;	
			flowData.capFlows(maxFlows);	
			
			Sequence currSeq = flowData.getSequence();
			
			if(!flowData.hasMinFlows(minFlows)){	//screen to see if sequence is of a minimum number of flows
				success = 0;
				trashCode += 'l';
			}
			
			int primerIndex = 0;
			int barcodeIndex = 0;
			
			if(barcodes.size() != 0){
				success = trimOligos.stripBarcode(currSeq, barcodeIndex);
				if(success > bdiffs)		{	trashCode += 'b';	}
				else{ currentSeqDiffs += success;  }
			}
			
			if(numFPrimers != 0){
				success = trimOligos.stripForward(currSeq, primerIndex);
				if(success > pdiffs)		{	trashCode += 'f';	}
				else{ currentSeqDiffs += success;  }
			}
			
			if (currentSeqDiffs > tdiffs)	{	trashCode += 't';   }
			
			if(numRPrimers != 0){
				success = trimOligos.stripReverse(currSeq);
				if(!success)				{	trashCode += 'r';	}
			}
			
			if(trashCode.length() == 0){
							
				flowData.printFlows(trimFlowFile);
			
				if(fasta)	{	currSeq.printSequence(fastaFile);	}
				
				if(allFiles){
					ofstream output;
					m->openOutputFileAppend(thisBarcodePrimerComboFileNames[barcodeIndex][primerIndex], output);
					output.setf(ios::fixed, ios::floatfield); trimFlowFile.setf(ios::showpoint);
					
					flowData.printFlows(output);
					output.close();
				}				
			}
			else{
				flowData.printFlows(scrapFlowFile, trashCode);
			}
				
			count++;
			//cout << "driver" << '\t' << currSeq.getName() << endl;			
			//report progress
			if((count) % 10000 == 0){	m->mothurOut(toString(count)); m->mothurOutEndLine();		}

#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
			unsigned long long pos = flowFile.tellg();

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
/**************************************************************************************************/
vector<unsigned long long> TrimFlowsCommand::getFlowFileBreaks() {

	try{
			
		vector<unsigned long long> filePos;
		filePos.push_back(0);
					
		FILE * pFile;
		unsigned long long size;
		
		//get num bytes in file
		pFile = fopen (flowFileName.c_str(),"rb");
		if (pFile==NULL) perror ("Error opening file");
		else{
			fseek (pFile, 0, SEEK_END);
			size=ftell (pFile);
			fclose (pFile);
		}
				
		//estimate file breaks
		unsigned long long chunkSize = 0;
		chunkSize = size / processors;

		//file too small to divide by processors
		if (chunkSize == 0)  {  processors = 1;	filePos.push_back(size); return filePos;	}
		
		//for each process seekg to closest file break and search for next '>' char. make that the filebreak
		for (int i = 0; i < processors; i++) {
			unsigned long long spot = (i+1) * chunkSize;
			
			ifstream in;
			m->openInputFile(flowFileName, in);
			in.seekg(spot);
			
			string dummy = m->getline(in);
			
			//there was not another sequence before the end of the file
			unsigned long long sanityPos = in.tellg();
			
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
		//unsigned long long spot = in.tellg();
		//filePos[0] = spot;
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
		processIDS.clear();
		int exitCommand = 1;
		
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
		int process = 1;
		
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
#else
		//////////////////////////////////////////////////////////////////////////////////////////////////////
		//Windows version shared memory, so be careful when passing variables through the trimFlowData struct. 
		//Above fork() will clone, so memory is separate, but that's not the case with windows, 
		//////////////////////////////////////////////////////////////////////////////////////////////////////
		
		vector<trimFlowData*> pDataArray; 
		DWORD   dwThreadIdArray[processors-1];
		HANDLE  hThreadArray[processors-1]; 
		
		//Create processor worker threads.
		for( int i=0; i<processors-1; i++ ){
			// Allocate memory for thread data.
			string extension = "";
			if (i != 0) { extension = toString(i) + ".temp"; processIDS.push_back(i); }
			
			vector<vector<string> > tempBarcodePrimerComboFileNames = barcodePrimerComboFileNames;
			if(allFiles){
				for(int i=0;i<tempBarcodePrimerComboFileNames.size();i++){
					for(int j=0;j<tempBarcodePrimerComboFileNames[0].size();j++){
						tempBarcodePrimerComboFileNames[i][j] += extension;
						ofstream temp;
						m->openOutputFile(tempBarcodePrimerComboFileNames[i][j], temp);
						temp.close();
						
					}
				}
			}
			
			trimFlowData* tempflow = new trimFlowData(flowFileName, (trimFlowFileName + extension), (scrapFlowFileName + extension), fastaFileName, flowOrder, tempBarcodePrimerComboFileNames, barcodes, primers, revPrimer, fasta, allFiles, lines[i]->start, lines[i]->end, m, signal, noise, numFlows, maxFlows, minFlows, maxHomoP, tdiffs, bdiffs, pdiffs, i);
			pDataArray.push_back(tempflow);
			
			//MyTrimFlowThreadFunction is in header. It must be global or static to work with the threads.
			//default security attributes, thread function name, argument to thread function, use default creation flags, returns the thread identifier
			hThreadArray[i] = CreateThread(NULL, 0, MyTrimFlowThreadFunction, pDataArray[i], 0, &dwThreadIdArray[i]);   
		}
		
		//using the main process as a worker saves time and memory
		ofstream temp;
		m->openOutputFile(trimFlowFileName, temp);
		temp.close();
		
		m->openOutputFile(scrapFlowFileName, temp);
		temp.close();
		
		if(fasta){
			m->openOutputFile(fastaFileName, temp);
			temp.close();
		}
		
		vector<vector<string> > tempBarcodePrimerComboFileNames = barcodePrimerComboFileNames;
		if(allFiles){
			for(int i=0;i<tempBarcodePrimerComboFileNames.size();i++){
				for(int j=0;j<tempBarcodePrimerComboFileNames[0].size();j++){
					tempBarcodePrimerComboFileNames[i][j] += toString(processors-1) + ".temp";
					ofstream temp;
					m->openOutputFile(tempBarcodePrimerComboFileNames[i][j], temp);
					temp.close();
					
				}
			}
		}
		
		//do my part - do last piece because windows is looking for eof
		int num = driverCreateTrim(flowFileName, (trimFlowFileName  + toString(processors-1) + ".temp"), (scrapFlowFileName  + toString(processors-1) + ".temp"), (fastaFileName + toString(processors-1) + ".temp"), tempBarcodePrimerComboFileNames, lines[processors-1]);
		processIDS.push_back((processors-1)); 
		
		//Wait until all threads have terminated.
		WaitForMultipleObjects(processors-1, hThreadArray, TRUE, INFINITE);
		
		//Close all thread handles and free memory allocations.
		for(int i=0; i < pDataArray.size(); i++){
			num += pDataArray[i]->count;
			CloseHandle(hThreadArray[i]);
			delete pDataArray[i];
		}
		
		
#endif	
		//append files
		m->mothurOutEndLine();
		for(int i=0;i<processIDS.size();i++){
			
			m->mothurOut("Appending files from process " + toString(processIDS[i])); m->mothurOutEndLine();
			
			m->appendFiles((trimFlowFileName + toString(processIDS[i]) + ".temp"), trimFlowFileName);
			m->mothurRemove((trimFlowFileName + toString(processIDS[i]) + ".temp"));
//			m->mothurOut("\tDone with trim.flow file"); m->mothurOutEndLine();

			m->appendFiles((scrapFlowFileName + toString(processIDS[i]) + ".temp"), scrapFlowFileName);
			m->mothurRemove((scrapFlowFileName + toString(processIDS[i]) + ".temp"));
//			m->mothurOut("\tDone with scrap.flow file"); m->mothurOutEndLine();

			if(fasta){
				m->appendFiles((fastaFileName + toString(processIDS[i]) + ".temp"), fastaFileName);
				m->mothurRemove((fastaFileName + toString(processIDS[i]) + ".temp"));
//				m->mothurOut("\tDone with flow.fasta file"); m->mothurOutEndLine();
			}
			if(allFiles){						
				for (int j = 0; j < barcodePrimerComboFileNames.size(); j++) {
					for (int k = 0; k < barcodePrimerComboFileNames[0].size(); k++) {
						m->appendFiles((barcodePrimerComboFileNames[j][k] + toString(processIDS[i]) + ".temp"), barcodePrimerComboFileNames[j][k]);
						m->mothurRemove((barcodePrimerComboFileNames[j][k] + toString(processIDS[i]) + ".temp"));
					}
				}
			}
		}
		
		return exitCommand;
	
	}
	catch(exception& e) {
		m->errorOut(e, "TrimFlowsCommand", "createProcessesCreateTrim");
		exit(1);
	}
}

//***************************************************************************************************************
