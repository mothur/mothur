/*
 *  trimseqscommand.cpp
 *  Mothur
 *
 *  Created by Pat Schloss on 6/6/09.
 *  Copyright 2009 Patrick D. Schloss. All rights reserved.
 *
 */

#include "trimseqscommand.h"
#include "needlemanoverlap.hpp"

//**********************************************************************************************************************

vector<string> TrimSeqsCommand::getValidParameters(){	
	try {
		string Array[] =  {"fasta", "flip", "oligos", "maxambig", "maxhomop","minlength", "maxlength", "qfile", 
									"qthreshold", "qwindowaverage", "qstepsize", "qwindowsize", "qaverage", "rollaverage",
									"keepfirst", "removelast",
									"allfiles", "qtrim","tdiffs", "pdiffs", "bdiffs", "processors", "outputdir","inputdir"};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "TrimSeqsCommand", "getValidParameters");
		exit(1);
	}
}

//**********************************************************************************************************************

TrimSeqsCommand::TrimSeqsCommand(){	
	try {
		abort = true; calledHelp = true; 
		vector<string> tempOutNames;
		outputTypes["fasta"] = tempOutNames;
		outputTypes["qfile"] = tempOutNames;
		outputTypes["group"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "TrimSeqsCommand", "TrimSeqsCommand");
		exit(1);
	}
}

//**********************************************************************************************************************

vector<string> TrimSeqsCommand::getRequiredParameters(){	
	try {
		string Array[] =  {"fasta"};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "TrimSeqsCommand", "getRequiredParameters");
		exit(1);
	}
}

//**********************************************************************************************************************

vector<string> TrimSeqsCommand::getRequiredFiles(){	
	try {
		vector<string> myArray;
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "TrimSeqsCommand", "getRequiredFiles");
		exit(1);
	}
}

//***************************************************************************************************************

TrimSeqsCommand::TrimSeqsCommand(string option)  {
	try {
		
		abort = false; calledHelp = false;   
		comboStarts = 0;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		
		else {
			//valid paramters for this command
			string AlignArray[] =  {	"fasta", "flip", "oligos", "maxambig", "maxhomop", "minlength", "maxlength", "qfile", 
								"qthreshold", "qwindowaverage", "qstepsize", "qwindowsize", "qaverage", "rollaverage",
								"keepfirst", "removelast",
								"allfiles", "qtrim","tdiffs", "pdiffs", "bdiffs", "processors", "outputdir","inputdir"};
			
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
			outputTypes["fasta"] = tempOutNames;
			outputTypes["qfile"] = tempOutNames;
			outputTypes["group"] = tempOutNames;
			
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
				
				it = parameters.find("oligos");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["oligos"] = inputDir + it->second;		}
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
			fastaFile = validParameter.validFile(parameters, "fasta", true);
			if (fastaFile == "not found") { m->mothurOut("fasta is a required parameter for the trim.seqs command."); m->mothurOutEndLine(); abort = true; }
			else if (fastaFile == "not open") { abort = true; }	
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += m->hasPath(fastaFile); //if user entered a file with a path then preserve it	
			}
		
			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			string temp;
			temp = validParameter.validFile(parameters, "flip", false);
			if (temp == "not found"){	flip = 0;	}
			else if(m->isTrue(temp))	{	flip = 1;	}
		
			temp = validParameter.validFile(parameters, "oligos", true);
			if (temp == "not found"){	oligoFile = "";		}
			else if(temp == "not open"){	abort = true;	} 
			else					{	oligoFile = temp;		}
			
			
			temp = validParameter.validFile(parameters, "maxambig", false);		if (temp == "not found") { temp = "-1"; }
			convert(temp, maxAmbig);  

			temp = validParameter.validFile(parameters, "maxhomop", false);		if (temp == "not found") { temp = "0"; }
			convert(temp, maxHomoP);  

			temp = validParameter.validFile(parameters, "minlength", false);	if (temp == "not found") { temp = "0"; }
			convert(temp, minLength); 
			
			temp = validParameter.validFile(parameters, "maxlength", false);	if (temp == "not found") { temp = "0"; }
			convert(temp, maxLength);
			
			temp = validParameter.validFile(parameters, "bdiffs", false);		if (temp == "not found") { temp = "0"; }
			convert(temp, bdiffs);
			
			temp = validParameter.validFile(parameters, "pdiffs", false);		if (temp == "not found") { temp = "0"; }
			convert(temp, pdiffs);
			
			temp = validParameter.validFile(parameters, "tdiffs", false);		if (temp == "not found") { int tempTotal = pdiffs + bdiffs;  temp = toString(tempTotal); }
			convert(temp, tdiffs);
			
			if(tdiffs == 0){	tdiffs = bdiffs + pdiffs;	}
			
			temp = validParameter.validFile(parameters, "qfile", true);	
			if (temp == "not found")	{	qFileName = "";		}
			else if(temp == "not open")	{	abort = true;		}
			else						{	qFileName = temp;	}
			
			temp = validParameter.validFile(parameters, "qthreshold", false);	if (temp == "not found") { temp = "0"; }
			convert(temp, qThreshold);
			
			temp = validParameter.validFile(parameters, "qtrim", false);		if (temp == "not found") { temp = "t"; }
			qtrim = m->isTrue(temp);

			temp = validParameter.validFile(parameters, "rollaverage", false);	if (temp == "not found") { temp = "0"; }
			convert(temp, qRollAverage);

			temp = validParameter.validFile(parameters, "qwindowaverage", false);if (temp == "not found") { temp = "0"; }
			convert(temp, qWindowAverage);

			temp = validParameter.validFile(parameters, "qwindowsize", false);	if (temp == "not found") { temp = "50"; }
			convert(temp, qWindowSize);

			temp = validParameter.validFile(parameters, "qstepsize", false);	if (temp == "not found") { temp = "1"; }
			convert(temp, qWindowStep);

			temp = validParameter.validFile(parameters, "qaverage", false);		if (temp == "not found") { temp = "0"; }
			convert(temp, qAverage);

			temp = validParameter.validFile(parameters, "keepfirst", false);	if (temp == "not found") { temp = "0"; }
			convert(temp, keepFirst);

			temp = validParameter.validFile(parameters, "removelast", false);	if (temp == "not found") { temp = "0"; }
			convert(temp, removeLast);
			
			temp = validParameter.validFile(parameters, "allfiles", false);		if (temp == "not found") { temp = "F"; }
			allFiles = m->isTrue(temp);
			
			temp = validParameter.validFile(parameters, "processors", false);	if (temp == "not found") { temp = "1"; }
			convert(temp, processors); 
			
			
			if(allFiles && (oligoFile == "")){
				m->mothurOut("You selected allfiles, but didn't enter an oligos.  Ignoring the allfiles request."); m->mothurOutEndLine();
			}
			if((qAverage != 0 && qThreshold != 0) && qFileName == ""){
				m->mothurOut("You didn't provide a quality file name, quality criteria will be ignored."); m->mothurOutEndLine();
				qAverage=0;
				qThreshold=0;
			}
			if(!flip && oligoFile=="" && !maxLength && !minLength && (maxAmbig==-1) && !maxHomoP && qFileName == ""){		
				m->mothurOut("You didn't set any options... quiting command."); m->mothurOutEndLine();
				abort = true;
			}
		}

	}
	catch(exception& e) {
		m->errorOut(e, "TrimSeqsCommand", "TrimSeqsCommand");
		exit(1);
	}
}

//**********************************************************************************************************************

void TrimSeqsCommand::help(){
	try {
		m->mothurOut("The trim.seqs command reads a fastaFile and creates 2 new fasta files, .trim.fasta and scrap.fasta, as well as group files if you provide and oligos file.\n");
		m->mothurOut("The .trim.fasta contains sequences that meet your requirements, and the .scrap.fasta contains those which don't.\n");
		m->mothurOut("The trim.seqs command parameters are fasta, flip, oligos, maxambig, maxhomop, minlength, maxlength, qfile, qthreshold, qaverage, diffs, qtrim, keepfirst, removelast and allfiles.\n");
		m->mothurOut("The fasta parameter is required.\n");
		m->mothurOut("The flip parameter will output the reverse compliment of your trimmed sequence. The default is false.\n");
		m->mothurOut("The oligos parameter allows you to provide an oligos file.\n");
		m->mothurOut("The maxambig parameter allows you to set the maximum number of ambigious bases allowed. The default is -1.\n");
		m->mothurOut("The maxhomop parameter allows you to set a maximum homopolymer length. \n");
		m->mothurOut("The minlength parameter allows you to set and minimum sequence length. \n");
		m->mothurOut("The maxlength parameter allows you to set and maximum sequence length. \n");
		m->mothurOut("The tdiffs parameter is used to specify the total number of differences allowed in the sequence. The default is pdiffs + bdiffs.\n");
		m->mothurOut("The bdiffs parameter is used to specify the number of differences allowed in the barcode. The default is 0.\n");
		m->mothurOut("The pdiffs parameter is used to specify the number of differences allowed in the primer. The default is 0.\n");
		m->mothurOut("The qfile parameter allows you to provide a quality file.\n");
		m->mothurOut("The qthreshold parameter allows you to set a minimum quality score allowed. \n");
		m->mothurOut("The qaverage parameter allows you to set a minimum average quality score allowed. \n");
		m->mothurOut("The qwindowsize parameter allows you to set a number of bases in a window. Default=50.\n");
		m->mothurOut("The qwindowaverage parameter allows you to set a minimum average quality score allowed over a window. \n");
		m->mothurOut("The rollaverage parameter allows you to set a minimum rolling average quality score allowed over a window. \n");
		m->mothurOut("The qstepsize parameter allows you to set a number of bases to move the window over. Default=1.\n");
		m->mothurOut("The allfiles parameter will create separate group and fasta file for each grouping. The default is F.\n");
		m->mothurOut("The qtrim parameter will trim sequence from the point that they fall below the qthreshold and put it in the .trim file if set to true. The default is T.\n");
		m->mothurOut("The keepfirst parameter trims the sequence to the first keepfirst number of bases after the barcode or primers are removed, before the sequence is checked to see if it meets the other requirements. \n");
		m->mothurOut("The removelast removes the last removelast number of bases after the barcode or primers are removed, before the sequence is checked to see if it meets the other requirements.\n");
		m->mothurOut("The trim.seqs command should be in the following format: \n");
		m->mothurOut("trim.seqs(fasta=yourFastaFile, flip=yourFlip, oligos=yourOligos, maxambig=yourMaxambig,  \n");
		m->mothurOut("maxhomop=yourMaxhomop, minlength=youMinlength, maxlength=yourMaxlength)  \n");	
		m->mothurOut("Example trim.seqs(fasta=abrecovery.fasta, flip=..., oligos=..., maxambig=..., maxhomop=..., minlength=..., maxlength=...).\n");
		m->mothurOut("Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFasta).\n");
		m->mothurOut("For more details please check out the wiki http://www.mothur.org/wiki/Trim.seqs .\n\n");

	}
	catch(exception& e) {
		m->errorOut(e, "TrimSeqsCommand", "help");
		exit(1);
	}
}


//***************************************************************************************************************

TrimSeqsCommand::~TrimSeqsCommand(){	/*	do nothing	*/	}

//***************************************************************************************************************

int TrimSeqsCommand::execute(){
	try{
	
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		numFPrimers = 0;  //this needs to be initialized
		numRPrimers = 0;
		vector<vector<string> > fastaFileNames;
		vector<vector<string> > qualFileNames;
		
		string trimSeqFile = outputDir + m->getRootName(m->getSimpleName(fastaFile)) + "trim.fasta";
		outputNames.push_back(trimSeqFile); outputTypes["fasta"].push_back(trimSeqFile);
		
		string scrapSeqFile = outputDir + m->getRootName(m->getSimpleName(fastaFile)) + "scrap.fasta";
		outputNames.push_back(scrapSeqFile); outputTypes["fasta"].push_back(scrapSeqFile);
		
		string trimQualFile = outputDir + m->getRootName(m->getSimpleName(fastaFile)) + "trim.qual";
		string scrapQualFile = outputDir + m->getRootName(m->getSimpleName(fastaFile)) + "scrap.qual";
		if (qFileName != "") {
			outputNames.push_back(trimQualFile);
			outputNames.push_back(scrapQualFile);
			outputTypes["qfile"].push_back(trimQualFile);
			outputTypes["qfile"].push_back(scrapQualFile); 
		}
		
		string outputGroupFileName;

		if(oligoFile != ""){
			outputGroupFileName = outputDir + m->getRootName(m->getSimpleName(fastaFile)) + "groups";
			outputNames.push_back(outputGroupFileName); outputTypes["group"].push_back(outputGroupFileName);
			getOligos(fastaFileNames, qualFileNames);
		}

		vector<unsigned long int> fastaFilePos;
		vector<unsigned long int> qFilePos;
		
		setLines(fastaFile, qFileName, fastaFilePos, qFilePos);
		
		for (int i = 0; i < (fastaFilePos.size()-1); i++) {
			lines.push_back(new linePair(fastaFilePos[i], fastaFilePos[(i+1)]));
			if (qFileName != "") {  qLines.push_back(new linePair(qFilePos[i], qFilePos[(i+1)]));  }
		}	
		if(qFileName == "")	{	qLines = lines;	} //files with duds
		
		#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
				if(processors == 1){
					driverCreateTrim(fastaFile, qFileName, trimSeqFile, scrapSeqFile, trimQualFile, scrapQualFile, outputGroupFileName, fastaFileNames, qualFileNames, lines[0], qLines[0]);
				}else{
					createProcessesCreateTrim(fastaFile, qFileName, trimSeqFile, scrapSeqFile, trimQualFile, scrapQualFile, outputGroupFileName, fastaFileNames, qualFileNames); 
				}	
		#else
				driverCreateTrim(fastaFile, qFileName, trimSeqFile, scrapSeqFile, trimQualFile, scrapQualFile, outputGroupFileName, fastaFileNames, qualFileNames, lines[0], qLines[0]);
		#endif
		
		if (m->control_pressed) {  return 0; }			
			
				
		if(allFiles){
			for(int i=0;i<fastaFileNames.size();i++){
				for(int j=0;j<fastaFileNames[0].size();j++){
					if(m->isBlank(fastaFileNames[i][j])){
						remove(fastaFileNames[i][j].c_str());

						if(qFileName != ""){
							remove(fastaFileNames[i][j].c_str());
						}

					}
				}
			}
		}
		
		
		
		if (m->control_pressed) { 
			for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); }
			return 0;
		}

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
		
/**************************************************************************************/

int TrimSeqsCommand::driverCreateTrim(string filename, string qFileName, string trimFileName, string scrapFileName, string trimQFileName, string scrapQFileName, string groupFileName, vector<vector<string> > fastaFileNames, vector<vector<string> > qualFileNames, linePair* line, linePair* qline) {	
		
	try {
		
		ofstream trimFASTAFile;
		m->openOutputFile(trimFileName, trimFASTAFile);
		
		ofstream scrapFASTAFile;
		m->openOutputFile(scrapFileName, scrapFASTAFile);
		
		ofstream trimQualFile;
		ofstream scrapQualFile;
		if(qFileName != ""){
			m->openOutputFile(trimQFileName, trimQualFile);
			m->openOutputFile(scrapQFileName, scrapQualFile);
		}
		
		ofstream outGroupsFile;
		if (oligoFile != ""){	m->openOutputFile(groupFileName, outGroupsFile);   }
		
		if(allFiles){
			for (int i = 0; i < fastaFileNames.size(); i++) { //clears old file
				for (int j = 0; j < fastaFileNames[i].size(); j++) { //clears old file
					ofstream temp;
					m->openOutputFile(fastaFileNames[i][j], temp);			temp.close();
					if(qFileName != ""){
						m->openOutputFile(qualFileNames[i][j], temp);			temp.close();
					}
				}
			}
		}
		
		ifstream inFASTA;
		m->openInputFile(filename, inFASTA);
		inFASTA.seekg(line->start);
		
		ifstream qFile;
		if(qFileName != "")	{
			m->openInputFile(qFileName, qFile);
			qFile.seekg(qline->start);  
		}
		
		int count = 0;
		bool moreSeqs = 1;
	
		while (moreSeqs) {
				
			if (m->control_pressed) { 
				inFASTA.close(); trimFASTAFile.close(); scrapFASTAFile.close();
				if (oligoFile != "") {	 outGroupsFile.close();   }

				if(qFileName != ""){
					qFile.close();
				}
				for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); }

				return 0;
			}
			
			int success = 1;
			string trashCode = "";
			int currentSeqsDiffs = 0;

			Sequence currSeq(inFASTA); m->gobble(inFASTA);

			QualityScores currQual;
			if(qFileName != ""){
				currQual = QualityScores(qFile);  m->gobble(qFile);
			}

			string origSeq = currSeq.getUnaligned();
			if (origSeq != "") {
				
				int barcodeIndex = 0;
				int primerIndex = 0;
				
				if(barcodes.size() != 0){
					success = stripBarcode(currSeq, currQual, barcodeIndex);
					if(success > bdiffs)		{	trashCode += 'b';	}
					else{ currentSeqsDiffs += success;  }
				}
				
				if(numFPrimers != 0){
					success = stripForward(currSeq, currQual, primerIndex);
					if(success > pdiffs)		{	trashCode += 'f';	}
					else{ currentSeqsDiffs += success;  }
				}
				
				if (currentSeqsDiffs > tdiffs)	{	trashCode += 't';   }
				
				if(numRPrimers != 0){
					success = stripReverse(currSeq, currQual);
					if(!success)				{	trashCode += 'r';	}
				}

				if(keepFirst != 0){
					success = keepFirstTrim(currSeq, currQual);
				}
				
				if(removeLast != 0){
					success = removeLastTrim(currSeq, currQual);
					if(!success)				{	trashCode += 'l';	}
				}

				
				if(qFileName != ""){
					int origLength = currSeq.getNumBases();
					
					if(qThreshold != 0)			{	success = currQual.stripQualThreshold(currSeq, qThreshold);			}
					else if(qAverage != 0)		{	success = currQual.cullQualAverage(currSeq, qAverage);				}
					else if(qRollAverage != 0)	{	success = currQual.stripQualRollingAverage(currSeq, qRollAverage);	}
					else if(qWindowAverage != 0){	success = currQual.stripQualWindowAverage(currSeq, qWindowStep, qWindowSize, qWindowAverage);	}
					else						{	success = 1;				}
					
					//you don't want to trim, if it fails above then scrap it
					if ((!qtrim) && (origLength != currSeq.getNumBases())) {  success = 0; }
					
					if(!success)				{	trashCode += 'q';	}
				}				
		
				if(minLength > 0 || maxLength > 0){
					success = cullLength(currSeq);
					if(!success)				{	trashCode += 'l';	}
				}
				if(maxHomoP > 0){
					success = cullHomoP(currSeq);
					if(!success)				{	trashCode += 'h';	}
				}
				if(maxAmbig != -1){
					success = cullAmbigs(currSeq);
					if(!success)				{	trashCode += 'n';	}
				}
				
				if(flip){		// should go last			
					currSeq.reverseComplement();
					if(qFileName != ""){
						currQual.flipQScores();	
					}
				}
				
				if(trashCode.length() == 0){
					currSeq.setAligned(currSeq.getUnaligned());
					currSeq.printSequence(trimFASTAFile);
					
					if(qFileName != ""){
						currQual.printQScores(trimQualFile);
					}
					
					if(barcodes.size() != 0){
						outGroupsFile << currSeq.getName() << '\t' << barcodeNameVector[barcodeIndex] << endl;
					}
					
					
					if(allFiles){
						ofstream output;
						m->openOutputFileAppend(fastaFileNames[barcodeIndex][primerIndex], output);
						currSeq.printSequence(output);
						output.close();
						
						if(qFileName != ""){
							m->openOutputFileAppend(qualFileNames[barcodeIndex][primerIndex], output);
							currQual.printQScores(output);
							output.close();							
						}
					}
				}
				else{
					currSeq.setName(currSeq.getName() + '|' + trashCode);
					currSeq.setUnaligned(origSeq);
					currSeq.setAligned(origSeq);
					currSeq.printSequence(scrapFASTAFile);
					if(qFileName != ""){
						currQual.printQScores(scrapQualFile);
					}
				}
				count++;
			}
			
			#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
				unsigned long int pos = inFASTA.tellg();
				if ((pos == -1) || (pos >= line->end)) { break; }
			#else
				if (inFASTA.eof()) { break; }
			#endif
				
			//report progress
			if((count) % 1000 == 0){	m->mothurOut(toString(count)); m->mothurOutEndLine();		}
			
		}
		//report progress
		if((count) % 1000 != 0){	m->mothurOut(toString(count)); m->mothurOutEndLine();		}

		
		inFASTA.close();
		trimFASTAFile.close();
		scrapFASTAFile.close();
		if (oligoFile != "") {	 outGroupsFile.close();   }
		if(qFileName != "")	{	qFile.close();	scrapQualFile.close(); trimQualFile.close();	}
		
		return count;
	}
	catch(exception& e) {
		m->errorOut(e, "TrimSeqsCommand", "driverCreateTrim");
		exit(1);
	}
}

/**************************************************************************************************/

int TrimSeqsCommand::createProcessesCreateTrim(string filename, string qFileName, string trimFASTAFileName, string scrapFASTAFileName, string trimQualFileName, string scrapQualFileName, string groupFile, vector<vector<string> > fastaFileNames, vector<vector<string> > qualFileNames) {
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
				
				vector<vector<string> > tempFASTAFileNames = fastaFileNames;
				vector<vector<string> > tempPrimerQualFileNames = qualFileNames;

				if(allFiles){
					ofstream temp;

					for(int i=0;i<tempFASTAFileNames.size();i++){
						for(int j=0;j<tempFASTAFileNames[i].size();j++){
							tempFASTAFileNames[i][j] += toString(getpid()) + ".temp";
							m->openOutputFile(tempFASTAFileNames[i][j], temp);			temp.close();

							if(qFileName != ""){
								tempPrimerQualFileNames[i][j] += toString(getpid()) + ".temp";
								m->openOutputFile(tempPrimerQualFileNames[i][j], temp);		temp.close();
							}
						}
					}
				}
							
				driverCreateTrim(filename,
								 qFileName,
								 (trimFASTAFileName + toString(getpid()) + ".temp"),
								 (scrapFASTAFileName + toString(getpid()) + ".temp"),
								 (trimQualFileName + toString(getpid()) + ".temp"),
								 (scrapQualFileName + toString(getpid()) + ".temp"),
								 (groupFile + toString(getpid()) + ".temp"),
								 tempFASTAFileNames,
								 tempPrimerQualFileNames,
								 lines[process],
								 qLines[process]);
				
				exit(0);
			}else { 
				m->mothurOut("[ERROR]: unable to spawn the necessary processes."); m->mothurOutEndLine(); 
				for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
				exit(0);
			}
		}
		
		//parent do my part
		ofstream temp;
		m->openOutputFile(trimFASTAFileName, temp);		temp.close();
		m->openOutputFile(scrapFASTAFileName, temp);	temp.close();
		m->openOutputFile(trimQualFileName, temp);		temp.close();
		m->openOutputFile(scrapQualFileName, temp);		temp.close();

		
		
		driverCreateTrim(filename, qFileName, trimFASTAFileName, scrapFASTAFileName, trimQualFileName, scrapQualFileName, groupFile, fastaFileNames, qualFileNames, lines[0], qLines[0]);
		
		
		//force parent to wait until all the processes are done
		for (int i=0;i<processIDS.size();i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
		
		//append files
		for(int i=0;i<processIDS.size();i++){
			
			m->mothurOut("Appending files from process " + toString(processIDS[i])); m->mothurOutEndLine();
			
			m->appendFiles((trimFASTAFileName + toString(processIDS[i]) + ".temp"), trimFASTAFileName);
			remove((trimFASTAFileName + toString(processIDS[i]) + ".temp").c_str());
			m->appendFiles((scrapFASTAFileName + toString(processIDS[i]) + ".temp"), scrapFASTAFileName);
			remove((scrapFASTAFileName + toString(processIDS[i]) + ".temp").c_str());
			
			if(qFileName != ""){
				m->appendFiles((trimQualFileName + toString(processIDS[i]) + ".temp"), trimQualFileName);
				remove((trimQualFileName + toString(processIDS[i]) + ".temp").c_str());
				m->appendFiles((scrapQualFileName + toString(processIDS[i]) + ".temp"), scrapQualFileName);
				remove((scrapQualFileName + toString(processIDS[i]) + ".temp").c_str());
			}
			
			m->appendFiles((groupFile + toString(processIDS[i]) + ".temp"), groupFile);
			remove((groupFile + toString(processIDS[i]) + ".temp").c_str());
			
			
			if(allFiles){
				for(int j=0;j<fastaFileNames.size();j++){
					for(int k=0;k<fastaFileNames[j].size();k++){
						m->appendFiles((fastaFileNames[j][k] + toString(processIDS[i]) + ".temp"), fastaFileNames[j][k]);
						remove((fastaFileNames[j][k] + toString(processIDS[i]) + ".temp").c_str());
						
						if(qFileName != ""){
							m->appendFiles((qualFileNames[j][k] + toString(processIDS[i]) + ".temp"), qualFileNames[j][k]);
							remove((qualFileNames[j][k] + toString(processIDS[i]) + ".temp").c_str());
						}
					}
				}
			}
			
		}
	
		return exitCommand;
#endif		
	}
	catch(exception& e) {
		m->errorOut(e, "TrimSeqsCommand", "createProcessesCreateTrim");
		exit(1);
	}
}

/**************************************************************************************************/

int TrimSeqsCommand::setLines(string filename, string qfilename, vector<unsigned long int>& fastaFilePos, vector<unsigned long int>& qfileFilePos) {
	try {
		
		//set file positions for fasta file
		fastaFilePos = m->divideFile(filename, processors);
		
		if (qfilename == "") { return processors; }
		
		//get name of first sequence in each chunk
		map<string, int> firstSeqNames;
		for (int i = 0; i < (fastaFilePos.size()-1); i++) {
			ifstream in;
			m->openInputFile(filename, in);
			in.seekg(fastaFilePos[i]);
		
			Sequence temp(in); 
			firstSeqNames[temp.getName()] = i;
		
			in.close();
		}
				
		//seach for filePos of each first name in the qfile and save in qfileFilePos
		ifstream inQual;
		m->openInputFile(qfilename, inQual);
		
		string input;
		while(!inQual.eof()){	
			input = m->getline(inQual);

			if (input.length() != 0) {
				if(input[0] == '>'){ //this is a sequence name line
					istringstream nameStream(input);
					
					string sname = "";  nameStream >> sname;
					sname = sname.substr(1);
					
					map<string, int>::iterator it = firstSeqNames.find(sname);
					
					if(it != firstSeqNames.end()) { //this is the start of a new chunk
						unsigned long int pos = inQual.tellg(); 
						qfileFilePos.push_back(pos - input.length() - 1);	
						firstSeqNames.erase(it);
					}
				}
			}
			
			if (firstSeqNames.size() == 0) { break; }
		}
		inQual.close();
		
		
		if (firstSeqNames.size() != 0) { 
			for (map<string, int>::iterator it = firstSeqNames.begin(); it != firstSeqNames.end(); it++) {
				m->mothurOut(it->first + " is in your fasta file and not in your quality file, not using quality file."); m->mothurOutEndLine();
			}
			qFileName = "";
			return processors;
		}

		//get last file position of qfile
		FILE * pFile;
		unsigned long int size;
		
		//get num bytes in file
		pFile = fopen (qfilename.c_str(),"rb");
		if (pFile==NULL) perror ("Error opening file");
		else{
			fseek (pFile, 0, SEEK_END);
			size=ftell (pFile);
			fclose (pFile);
		}
		
		qfileFilePos.push_back(size);
		
		return processors;
	}
	catch(exception& e) {
		m->errorOut(e, "TrimSeqsCommand", "setLines");
		exit(1);
	}
}

//***************************************************************************************************************

void TrimSeqsCommand::getOligos(vector<vector<string> >& fastaFileNames, vector<vector<string> >& qualFileNames){
	try {
		ifstream inOligos;
		m->openInputFile(oligoFile, inOligos);
		
		ofstream test;
		
		string type, oligo, group;

		int indexPrimer = 0;
		int indexBarcode = 0;
		
		while(!inOligos.eof()){

			inOligos >> type; m->gobble(inOligos);
		 			
			if(type[0] == '#'){
				while (!inOligos.eof())	{	char c = inOligos.get(); if (c == 10 || c == 13){	break;	}	} // get rest of line if there's any crap there
			}
			else{
				//make type case insensitive
				for(int i=0;i<type.length();i++){	type[i] = toupper(type[i]);  }
				
				inOligos >> oligo;
				
				for(int i=0;i<oligo.length();i++){
					oligo[i] = toupper(oligo[i]);
					if(oligo[i] == 'U')	{	oligo[i] = 'T';	}
				}
				
				if(type == "FORWARD"){
					group = "";
					
					// get rest of line in case there is a primer name
					while (!inOligos.eof())	{	
						char c = inOligos.get(); 
						if (c == 10 || c == 13){	break;	}
						else if (c == 32 || c == 9){;} //space or tab
						else { 	group += c;  }
					} 
					
					//check for repeat barcodes
					map<string, int>::iterator itPrime = primers.find(oligo);
					if (itPrime != primers.end()) { m->mothurOut("primer " + oligo + " is in your oligos file already."); m->mothurOutEndLine();  }
					
					primers[oligo]=indexPrimer; indexPrimer++;		
					primerNameVector.push_back(group);
				}
				else if(type == "REVERSE"){
					Sequence oligoRC("reverse", oligo);
					oligoRC.reverseComplement();
					revPrimer.push_back(oligoRC.getUnaligned());
				}
				else if(type == "BARCODE"){
					inOligos >> group;
					
					//check for repeat barcodes
					map<string, int>::iterator itBar = barcodes.find(oligo);
					if (itBar != barcodes.end()) { m->mothurOut("barcode " + oligo + " is in your oligos file already."); m->mothurOutEndLine();  }
						
					barcodes[oligo]=indexBarcode; indexBarcode++;
					barcodeNameVector.push_back(group);
				}
				else{	m->mothurOut(type + " is not recognized as a valid type. Choices are forward, reverse, and barcode. Ignoring " + oligo + "."); m->mothurOutEndLine();  }
			}
			m->gobble(inOligos);
		}	
		inOligos.close();
		
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
		
		fastaFileNames.resize(barcodeNameVector.size());
		for(int i=0;i<fastaFileNames.size();i++){
			fastaFileNames[i].assign(primerNameVector.size(), "");
		}
		if(qFileName != ""){	qualFileNames = fastaFileNames;	}
		
		if(allFiles){
			for(map<string, int>::iterator itBar = barcodes.begin();itBar != barcodes.end();itBar++){
				for(map<string, int>::iterator itPrimer = primers.begin();itPrimer != primers.end(); itPrimer++){
					
					string primerName = primerNameVector[itPrimer->second];
					string barcodeName = barcodeNameVector[itBar->second];
					
					string comboGroupName = "";
					string fastaFileName = "";
					string qualFileName = "";
					
					if(primerName == ""){
						comboGroupName = barcodeNameVector[itBar->second];
					}
					else{
						if(barcodeName == ""){
							comboGroupName = primerNameVector[itPrimer->second];
						}
						else{
							comboGroupName = barcodeNameVector[itBar->second] + "." + primerNameVector[itPrimer->second];
						}
					}

					ofstream temp;
					fastaFileName = outputDir + m->getRootName(m->getSimpleName(fastaFile)) + comboGroupName + ".fasta";
					outputNames.push_back(fastaFileName);
					outputTypes["fasta"].push_back(fastaFileName);
					fastaFileNames[itBar->second][itPrimer->second] = fastaFileName;
					m->openOutputFile(fastaFileName, temp);		temp.close();

					if(qFileName != ""){
						qualFileName = outputDir + m->getRootName(m->getSimpleName(qFileName)) + comboGroupName + ".qual";
						outputNames.push_back(qualFileName);
						outputTypes["qfile"].push_back(qualFileName);
						qualFileNames[itBar->second][itPrimer->second] = qualFileName;
						m->openOutputFile(qualFileName, temp);		temp.close();
					}
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

int TrimSeqsCommand::stripBarcode(Sequence& seq, QualityScores& qual, int& group){
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
				
				if(qual.getName() != ""){
					qual.trimQScores(oligo.length(), -1);
				}
				
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
				
				if(qual.getName() != ""){
					qual.trimQScores(minPos, -1);
				}
				success = minDiff;
			}
			
			if (alignment != NULL) {  delete alignment;  }
			
		}
		
		return success;
		
	}
	catch(exception& e) {
		m->errorOut(e, "TrimSeqsCommand", "stripBarcode");
		exit(1);
	}

}

//***************************************************************************************************************

int TrimSeqsCommand::stripForward(Sequence& seq, QualityScores& qual, int& group){
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
				if(qual.getName() != ""){
					qual.trimQScores(oligo.length(), -1);
				}
				success = 0;
				break;
			}
		}

		//if you found the barcode or if you don't want to allow for diffs
		if ((pdiffs == 0) || (success == 0)) { return success;  }
		
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
				if(qual.getName() != ""){
					qual.trimQScores(minPos, -1);
				}
				success = minDiff;
			}
			
			if (alignment != NULL) {  delete alignment;  }
			
		}
		
		return success;

	}
	catch(exception& e) {
		m->errorOut(e, "TrimSeqsCommand", "stripForward");
		exit(1);
	}
}

//***************************************************************************************************************

bool TrimSeqsCommand::stripReverse(Sequence& seq, QualityScores& qual){
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
				if(qual.getName() != ""){
					qual.trimQScores(-1, rawSequence.length()-oligo.length());
				}
				success = 1;
				break;
			}
		}	
		return success;
		
	}
	catch(exception& e) {
		m->errorOut(e, "TrimSeqsCommand", "stripReverse");
		exit(1);
	}
}

//***************************************************************************************************************

bool TrimSeqsCommand::keepFirstTrim(Sequence& sequence, QualityScores& qscores){
	try {
		bool success = 1;
		if(qscores.getName() != ""){
			qscores.trimQScores(-1, keepFirst);
		}
		sequence.trim(keepFirst);
		return success;
	}
	catch(exception& e) {
		m->errorOut(e, "keepFirstTrim", "countDiffs");
		exit(1);
	}
	
}	

//***************************************************************************************************************

bool TrimSeqsCommand::removeLastTrim(Sequence& sequence, QualityScores& qscores){
	try {
		bool success = 0;
		
		int length = sequence.getNumBases() - removeLast;
		
		if(length > 0){
			if(qscores.getName() != ""){
				qscores.trimQScores(-1, length);
			}
			sequence.trim(length);
			success = 1;
		}
		else{
			success = 0;
		}

		return success;
	}
	catch(exception& e) {
		m->errorOut(e, "removeLastTrim", "countDiffs");
		exit(1);
	}
	
}	

//***************************************************************************************************************

bool TrimSeqsCommand::cullLength(Sequence& seq){
	try {
	
		int length = seq.getNumBases();
		bool success = 0;	//guilty until proven innocent
		
		if(length >= minLength && maxLength == 0)			{	success = 1;	}
		else if(length >= minLength && length <= maxLength)	{	success = 1;	}
		else												{	success = 0;	}
		
		return success;
	
	}
	catch(exception& e) {
		m->errorOut(e, "TrimSeqsCommand", "cullLength");
		exit(1);
	}
	
}

//***************************************************************************************************************

bool TrimSeqsCommand::cullHomoP(Sequence& seq){
	try {
		int longHomoP = seq.getLongHomoPolymer();
		bool success = 0;	//guilty until proven innocent
		
		if(longHomoP <= maxHomoP){	success = 1;	}
		else					{	success = 0;	}
		
		return success;
	}
	catch(exception& e) {
		m->errorOut(e, "TrimSeqsCommand", "cullHomoP");
		exit(1);
	}
	
}

//***************************************************************************************************************

bool TrimSeqsCommand::cullAmbigs(Sequence& seq){
	try {
		int numNs = seq.getAmbigBases();
		bool success = 0;	//guilty until proven innocent
		
		if(numNs <= maxAmbig)	{	success = 1;	}
		else					{	success = 0;	}
		
		return success;
	}
	catch(exception& e) {
		m->errorOut(e, "TrimSeqsCommand", "cullAmbigs");
		exit(1);
	}
	
}

//***************************************************************************************************************

bool TrimSeqsCommand::compareDNASeq(string oligo, string seq){
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
		m->errorOut(e, "TrimSeqsCommand", "compareDNASeq");
		exit(1);
	}

}

//***************************************************************************************************************

int TrimSeqsCommand::countDiffs(string oligo, string seq){
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
		m->errorOut(e, "TrimSeqsCommand", "countDiffs");
		exit(1);
	}

}

//***************************************************************************************************************
