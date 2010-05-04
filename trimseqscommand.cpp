/*
 *  trimseqscommand.cpp
 *  Mothur
 *
 *  Created by Pat Schloss on 6/6/09.
 *  Copyright 2009 Patrick D. Schloss. All rights reserved.
 *
 */

#include "trimseqscommand.h"

//***************************************************************************************************************

TrimSeqsCommand::TrimSeqsCommand(string option)  {
	try {
		
		abort = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string AlignArray[] =  {"fasta", "flip", "oligos", "maxambig", "maxhomop", "minlength", "maxlength", "qfile", 
									"qthreshold", "qaverage", "allfiles", "qtrim", "processors", "outputdir","inputdir"};
			
			vector<string> myArray (AlignArray, AlignArray+(sizeof(AlignArray)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			map<string,string>::iterator it;
			
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
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
				
				it = parameters.find("oligos");
				//user has given a template file
				if(it != parameters.end()){ 
					path = hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["oligos"] = inputDir + it->second;		}
				}
				
				it = parameters.find("qfile");
				//user has given a template file
				if(it != parameters.end()){ 
					path = hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["qfile"] = inputDir + it->second;		}
				}
			}

			
			//check for required parameters
			fastaFile = validParameter.validFile(parameters, "fasta", true);
			if (fastaFile == "not found") { m->mothurOut("fasta is a required parameter for the screen.seqs command."); m->mothurOutEndLine(); abort = true; }
			else if (fastaFile == "not open") { abort = true; }	
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += hasPath(fastaFile); //if user entered a file with a path then preserve it	
			}
		
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			string temp;
			temp = validParameter.validFile(parameters, "flip", false);
			if (temp == "not found"){	flip = 0;	}
			else if(isTrue(temp))	{	flip = 1;	}
		
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
			
			temp = validParameter.validFile(parameters, "qfile", true);	
			if (temp == "not found")	{	qFileName = "";		}
			else if(temp == "not open")	{	abort = true;		}
			else						{	qFileName = temp;	}
			
			temp = validParameter.validFile(parameters, "qthreshold", false);	if (temp == "not found") { temp = "0"; }
			convert(temp, qThreshold);
			
			temp = validParameter.validFile(parameters, "qtrim", false);	if (temp == "not found") { temp = "F"; }
			qtrim = isTrue(temp);

			temp = validParameter.validFile(parameters, "qaverage", false);		if (temp == "not found") { temp = "0"; }
			convert(temp, qAverage);
			
			temp = validParameter.validFile(parameters, "allfiles", false);		if (temp == "not found") { temp = "F"; }
			allFiles = isTrue(temp);
			
			temp = validParameter.validFile(parameters, "processors", false);	if (temp == "not found"){	temp = "1";				}
			convert(temp, processors); 
			
			if(allFiles && oligoFile == ""){
				m->mothurOut("You selected allfiles, but didn't enter an oligos file.  Ignoring the allfiles request."); m->mothurOutEndLine();
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
		m->mothurOut("The trim.seqs command reads a fastaFile and creates .....\n");
		m->mothurOut("The trim.seqs command parameters are fasta, flip, oligos, maxambig, maxhomop, minlength, maxlength, qfile, qthreshold, qaverage, qtrim and allfiles.\n");
		m->mothurOut("The fasta parameter is required.\n");
		m->mothurOut("The flip parameter .... The default is 0.\n");
		m->mothurOut("The oligos parameter .... The default is "".\n");
		m->mothurOut("The maxambig parameter .... The default is -1.\n");
		m->mothurOut("The maxhomop parameter .... The default is 0.\n");
		m->mothurOut("The minlength parameter .... The default is 0.\n");
		m->mothurOut("The maxlength parameter .... The default is 0.\n");
		m->mothurOut("The qfile parameter .....\n");
		m->mothurOut("The qthreshold parameter .... The default is 0.\n");
		m->mothurOut("The qaverage parameter .... The default is 0.\n");
		m->mothurOut("The allfiles parameter .... The default is F.\n");
		m->mothurOut("The qtrim parameter .... The default is F.\n");
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
	
		if (abort == true) { return 0; }
		
		numFPrimers = 0;  //this needs to be initialized
		numRPrimers = 0;
		
		string trimSeqFile = outputDir + getRootName(getSimpleName(fastaFile)) + "trim.fasta";
		outputNames.push_back(trimSeqFile);
		string scrapSeqFile = outputDir + getRootName(getSimpleName(fastaFile)) + "scrap.fasta";
		outputNames.push_back(scrapSeqFile);
		string groupFile = outputDir + getRootName(getSimpleName(fastaFile)) + "groups";
		
		vector<string> fastaFileNames;
		if(oligoFile != ""){
			outputNames.push_back(groupFile);
			getOligos(fastaFileNames);
		}
		
		if(qFileName != "")	{	setLines(qFileName, qLines);	}

		
		#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
				if(processors == 1){
					ifstream inFASTA;
					openInputFile(fastaFile, inFASTA);
					int numSeqs=count(istreambuf_iterator<char>(inFASTA),istreambuf_iterator<char>(), '>');
					inFASTA.close();
					
					lines.push_back(new linePair(0, numSeqs));
					
					driverCreateTrim(fastaFile, qFileName, trimSeqFile, scrapSeqFile, groupFile, fastaFileNames, lines[0], lines[0]);
					
					for (int j = 0; j < fastaFileNames.size(); j++) {
						rename((fastaFileNames[j] + toString(getpid()) + ".temp").c_str(), fastaFileNames[j].c_str());
					}

				}else{
					setLines(fastaFile, lines);	
					if(qFileName == "")	{	qLines = lines;	}	
								
					createProcessesCreateTrim(fastaFile, qFileName, trimSeqFile, scrapSeqFile, groupFile, fastaFileNames); 
					
					rename((trimSeqFile + toString(processIDS[0]) + ".temp").c_str(), trimSeqFile.c_str());
					rename((scrapSeqFile + toString(processIDS[0]) + ".temp").c_str(), scrapSeqFile.c_str());
					rename((groupFile + toString(processIDS[0]) + ".temp").c_str(), groupFile.c_str());
					for (int j = 0; j < fastaFileNames.size(); j++) {
						rename((fastaFileNames[j] + toString(processIDS[0]) + ".temp").c_str(), fastaFileNames[j].c_str());
					}
					//append files
					for(int i=1;i<processors;i++){
						appendFiles((trimSeqFile + toString(processIDS[i]) + ".temp"), trimSeqFile);
						remove((trimSeqFile + toString(processIDS[i]) + ".temp").c_str());
						appendFiles((scrapSeqFile + toString(processIDS[i]) + ".temp"), scrapSeqFile);
						remove((scrapSeqFile + toString(processIDS[i]) + ".temp").c_str());
						appendFiles((groupFile + toString(processIDS[i]) + ".temp"), groupFile);
						remove((groupFile + toString(processIDS[i]) + ".temp").c_str());
						for (int j = 0; j < fastaFileNames.size(); j++) {
							appendFiles((fastaFileNames[j] + toString(processIDS[i]) + ".temp"), fastaFileNames[j]);
							remove((fastaFileNames[j] + toString(processIDS[i]) + ".temp").c_str());
						}
					}
				}
				
				if (m->control_pressed) {  return 0; }
		#else
				ifstream inFASTA;
				openInputFile(fastafileNames[s], inFASTA);
				numSeqs=count(istreambuf_iterator<char>(inFASTA),istreambuf_iterator<char>(), '>');
				inFASTA.close();
				
				lines.push_back(new linePair(0, numSeqs));
				
				driverCreateSummary(fastafile, qFileName, trimSeqFile, scrapSeqFile, groupFile, fastaFileNames, lines[0], lines[0]);
				
				if (m->control_pressed) {  return 0; }
		#endif
						
										
		for(int i=0;i<fastaFileNames.size();i++){
			ifstream inFASTA;
			string seqName;
			openInputFile(getRootName(fastaFile) + groupVector[i] + ".fasta", inFASTA);
			ofstream outGroups;
			openOutputFile(outputDir + getRootName(getSimpleName(fastaFile)) + groupVector[i] + ".groups", outGroups);
			outputNames.push_back(outputDir + getRootName(getSimpleName(fastaFile)) + groupVector[i] + ".groups");
			
			while(!inFASTA.eof()){
				if(inFASTA.get() == '>'){
					inFASTA >> seqName;
					outGroups << seqName << '\t' << groupVector[i] << endl;
				}
				while (!inFASTA.eof())	{	char c = inFASTA.get(); if (c == 10 || c == 13){	break;	}	}
			}
			outGroups.close();
			inFASTA.close();
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
int TrimSeqsCommand::driverCreateTrim(string filename, string qFileName, string trimFile, string scrapFile, string groupFile, vector<string> fastaNames, linePair* line, linePair* qline) {	
	try {
		
		ofstream outFASTA;
		int able = openOutputFile(trimFile, outFASTA);
		
		ofstream scrapFASTA;
		openOutputFile(scrapFile, scrapFASTA);
		
		ofstream outGroups;
		vector<ofstream*> fastaFileNames;
		if (oligoFile != "") {		
			openOutputFile(groupFile, outGroups);   
			for (int i = 0; i < fastaNames.size(); i++) {
				fastaFileNames.push_back(new ofstream((fastaNames[i] + toString(getpid()) + ".temp").c_str(), ios::ate)); 
			}
		}
		
		ifstream inFASTA;
		openInputFile(filename, inFASTA);
		
		ifstream qFile;
		if(qFileName != "")	{	openInputFile(qFileName, qFile);	}
		
		qFile.seekg(qline->start);
		inFASTA.seekg(line->start);
		
		for(int i=0;i<line->num;i++){
				
			if (m->control_pressed) { 
				inFASTA.close(); 
				outFASTA.close();
				scrapFASTA.close();
				if (oligoFile != "") {	 outGroups.close();   }
				if(qFileName != "")	{	qFile.close();	}
				for(int i=0;i<fastaFileNames.size();i++){
					fastaFileNames[i]->close();
					delete fastaFileNames[i];
				}	
				for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); }
				return 0;
			}
			
			bool success = 1;
			
			Sequence currSeq(inFASTA);

			string origSeq = currSeq.getUnaligned();
			if (origSeq != "") {
				int group;
				string trashCode = "";
				
				if(qFileName != ""){
					if(qThreshold != 0)		{	success = stripQualThreshold(currSeq, qFile);	}
					else if(qAverage != 0)	{	success = cullQualAverage(currSeq, qFile);		}
					if ((!qtrim) && (origSeq.length() != currSeq.getUnaligned().length())) { 
						success = 0; //if you don't want to trim and the sequence does not meet quality requirements, move to scrap
					}
					if(!success)			{	trashCode += 'q';								}
				}
			
				if(barcodes.size() != 0){
					success = stripBarcode(currSeq, group);
					if(!success){	trashCode += 'b';	}
				}
			
				if(numFPrimers != 0){
					success = stripForward(currSeq);
					if(!success){	trashCode += 'f';	}
				}
					
				if(numRPrimers != 0){
					success = stripReverse(currSeq);
					if(!success){	trashCode += 'r';	}
				}
		
				if(minLength > 0 || maxLength > 0){
					success = cullLength(currSeq);
					if(!success){	trashCode += 'l'; }
				}
				if(maxHomoP > 0){
					success = cullHomoP(currSeq);
					if(!success){	trashCode += 'h';	}
				}
				if(maxAmbig != -1){
					success = cullAmbigs(currSeq);
					if(!success){	trashCode += 'n';	}
				}
				
				if(flip){	currSeq.reverseComplement();	}		// should go last			
				
				if(trashCode.length() == 0){
					currSeq.setAligned(currSeq.getUnaligned());  //this is because of a modification we made to the sequence class to fix a bug.  all seqs have an aligned version, which is the version that gets printed.
					currSeq.printSequence(outFASTA);
					if(barcodes.size() != 0){
						outGroups << currSeq.getName() << '\t' << groupVector[group] << endl;
						
						if(allFiles){
							currSeq.printSequence(*fastaFileNames[group]);					
						}
					}
				}
				else{
					currSeq.setName(currSeq.getName() + '|' + trashCode);
					currSeq.setUnaligned(origSeq);
					currSeq.printSequence(scrapFASTA);
				}
			}
			gobble(inFASTA);
		}
		
		inFASTA.close();
		outFASTA.close();
		scrapFASTA.close();
		if (oligoFile != "") {	 outGroups.close();   }
		if(qFileName != "")	{	qFile.close();	}
		
		for(int i=0;i<fastaFileNames.size();i++){
			fastaFileNames[i]->close();
			delete fastaFileNames[i];
		}		
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "TrimSeqsCommand", "driverCreateTrim");
		exit(1);
	}
}
/**************************************************************************************************/
int TrimSeqsCommand::createProcessesCreateTrim(string filename, string qFileName, string trimFile, string scrapFile, string groupFile, vector<string> fastaNames) {
	try {
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
		int process = 0;
		int exitCommand = 1;
		processIDS.clear();
		
		//loop through and create all the processes you want
		while (process != processors) {
			int pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
				process++;
			}else if (pid == 0){
				driverCreateTrim(filename, qFileName, (trimFile + toString(getpid()) + ".temp"), (scrapFile + toString(getpid()) + ".temp"), (groupFile + toString(getpid()) + ".temp"), fastaNames, lines[process], qLines[process]);
				exit(0);
			}else { m->mothurOut("unable to spawn the necessary processes."); m->mothurOutEndLine(); exit(0); }
		}
		
		//force parent to wait until all the processes are done
		for (int i=0;i<processors;i++) { 
			int temp = processIDS[i];
			wait(&temp);
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

int TrimSeqsCommand::setLines(string filename, vector<linePair*>& lines) {
	try {
		
		lines.clear();
		
		vector<long int> positions;
		
		ifstream inFASTA;
		openInputFile(filename, inFASTA);
			
		string input;
		while(!inFASTA.eof()){	
			input = getline(inFASTA);

			if (input.length() != 0) {
				if(input[0] == '>'){ long int pos = inFASTA.tellg(); positions.push_back(pos - input.length() - 1);	}
			}
		}
		inFASTA.close();
		
		int numFastaSeqs = positions.size();
	
		FILE * pFile;
		long size;
		
		//get num bytes in file
		pFile = fopen (filename.c_str(),"rb");
		if (pFile==NULL) perror ("Error opening file");
		else{
			fseek (pFile, 0, SEEK_END);
			size=ftell (pFile);
			fclose (pFile);
		}
		
		int numSeqsPerProcessor = numFastaSeqs / processors;
		
		for (int i = 0; i < processors; i++) {

			long int startPos = positions[ i * numSeqsPerProcessor ];
			if(i == processors - 1){
				numSeqsPerProcessor = numFastaSeqs - i * numSeqsPerProcessor;
			}else{  
				long int myEnd = positions[ (i+1) * numSeqsPerProcessor ];
			}
			lines.push_back(new linePair(startPos, numSeqsPerProcessor));
		}
		
		return numFastaSeqs;
	}
	catch(exception& e) {
		m->errorOut(e, "TrimSeqsCommand", "setLines");
		exit(1);
	}
}
//***************************************************************************************************************

void TrimSeqsCommand::getOligos(vector<string>& outFASTAVec){ //vector<ofstream*>& outFASTAVec
	try {
		ifstream inOligos;
		openInputFile(oligoFile, inOligos);
		
		ofstream test;
		
		string type, oligo, group;
		int index=0;
		
		while(!inOligos.eof()){
			inOligos >> type;
			
			if(type[0] == '#'){
				while (!inOligos.eof())	{	char c = inOligos.get(); if (c == 10 || c == 13){	break;	}	} // get rest of line if there's any crap there
			}
			else{
				inOligos >> oligo;
				
				for(int i=0;i<oligo.length();i++){
					oligo[i] = toupper(oligo[i]);
					if(oligo[i] == 'U')	{	oligo[i] = 'T';	}
				}
				
				if(type == "forward"){
					forPrimer.push_back(oligo);
				}
				else if(type == "reverse"){
					Sequence oligoRC("reverse", oligo);
					oligoRC.reverseComplement();
					revPrimer.push_back(oligoRC.getUnaligned());
				}
				else if(type == "barcode"){
					inOligos >> group;
					barcodes[oligo]=index++;
					groupVector.push_back(group);
					
					if(allFiles){
						//outFASTAVec.push_back(new ofstream((outputDir + getRootName(getSimpleName(fastaFile)) + group + ".fasta").c_str(), ios::ate));
						outputNames.push_back((outputDir + getRootName(getSimpleName(fastaFile)) + group + ".fasta"));
						outFASTAVec.push_back((outputDir + getRootName(getSimpleName(fastaFile)) + group + ".fasta"));
					}
				}
			}
		}
		
		inOligos.close();
		
		numFPrimers = forPrimer.size();
		numRPrimers = revPrimer.size();
		
	}
	catch(exception& e) {
		m->errorOut(e, "TrimSeqsCommand", "getOligos");
		exit(1);
	}

}

//***************************************************************************************************************

bool TrimSeqsCommand::stripBarcode(Sequence& seq, int& group){
	try {
		string rawSequence = seq.getUnaligned();
		bool success = 0;	//guilty until proven innocent
		
		for(map<string,int>::iterator it=barcodes.begin();it!=barcodes.end();it++){
			string oligo = it->first;
			if(rawSequence.length() < oligo.length()){	//let's just assume that the barcodes are the same length
				success = 0;
				break;
			}
			
			if(compareDNASeq(oligo, rawSequence.substr(0,oligo.length()))){
				group = it->second;
				seq.setUnaligned(rawSequence.substr(oligo.length()));
				success = 1;
				break;
			}
		}
		return success;
		
	}
	catch(exception& e) {
		m->errorOut(e, "TrimSeqsCommand", "stripBarcode");
		exit(1);
	}

}

//***************************************************************************************************************

bool TrimSeqsCommand::stripForward(Sequence& seq){
	try {
		string rawSequence = seq.getUnaligned();
		bool success = 0;	//guilty until proven innocent
		
		for(int i=0;i<numFPrimers;i++){
			string oligo = forPrimer[i];
			
			if(rawSequence.length() < oligo.length()){
				success = 0;
				break;
			}
			
			if(compareDNASeq(oligo, rawSequence.substr(0,oligo.length()))){
				seq.setUnaligned(rawSequence.substr(oligo.length()));
				success = 1;
				break;
			}
		}
		
		return success;
		
	}
	catch(exception& e) {
		m->errorOut(e, "TrimSeqsCommand", "stripForward");
		exit(1);
	}
}

//***************************************************************************************************************

bool TrimSeqsCommand::stripReverse(Sequence& seq){
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
		m->errorOut(e, "TrimSeqsCommand", "stripReverse");
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
				if(oligo[i] == 'A' || oligo[i] == 'T' || oligo[i] == 'G' || oligo[i] == 'C')	{	success = 0;	}
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
				
				if(success == 0)	{	break;	}
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

bool TrimSeqsCommand::stripQualThreshold(Sequence& seq, ifstream& qFile){
	try {
		string rawSequence = seq.getUnaligned();
		int seqLength;  // = rawSequence.length();
		string name, temp, temp2;
		
		qFile >> name >> temp;
	
		splitAtEquals(temp2, temp); //separates length=242, temp=length, temp2=242
		convert(temp, seqLength); //converts string to int
	
		if (name.length() != 0) {  if(name.substr(1) != seq.getName())	{	m->mothurOut("sequence name mismatch btwn fasta and qual file"); m->mothurOutEndLine();	}  } 
		while (!qFile.eof())	{	char c = qFile.get(); if (c == 10 || c == 13){	break;	}	}
		
		int score;
		int end = seqLength;
		
		for(int i=0;i<seqLength;i++){
			qFile >> score;
			
			if(score <= qThreshold){
				end = i;
				break;
			}
		}
		for(int i=end+1;i<seqLength;i++){
			qFile >> score;
		}
		
		seq.setUnaligned(rawSequence.substr(0,end));
		
		return 1;
	}
	catch(exception& e) {
		m->errorOut(e, "TrimSeqsCommand", "stripQualThreshold");
		exit(1);
	}
}

//***************************************************************************************************************

bool TrimSeqsCommand::cullQualAverage(Sequence& seq, ifstream& qFile){
	try {
		string rawSequence = seq.getUnaligned();
		int seqLength = seq.getNumBases();
		bool success = 0;	//guilty until proven innocent
		string name;
		
		qFile >> name;
		if (name[0] == '>') {  if(name.substr(1) != seq.getName())	{	m->mothurOut("sequence name mismatch btwn fasta: " + seq.getName() + " and qual file: " + name); m->mothurOutEndLine();	}  }
		
		while (!qFile.eof())	{	char c = qFile.get(); if (c == 10 || c == 13){	break;	}	}
		
		float score;	
		float average = 0;
		
		for(int i=0;i<seqLength;i++){
			qFile >> score;
			average += score;
		}
		average /= seqLength;

		if(average >= qAverage)	{	success = 1;	}
		else					{	success = 0;	}
		
		return success;
	}
	catch(exception& e) {
		m->errorOut(e, "TrimSeqsCommand", "cullQualAverage");
		exit(1);
	}
}

//***************************************************************************************************************
