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

TrimSeqsCommand::TrimSeqsCommand(string option){
	try {
		
		abort = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string AlignArray[] =  {"fasta", "flip", "oligos", "maxambig", "maxhomop", "minlength", "maxlength", "qfile", "qthreshold", "qaverage", "allfiles"};
			
			vector<string> myArray (AlignArray, AlignArray+(sizeof(AlignArray)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			
			//check to make sure all parameters are valid for command
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//check for required parameters
			fastaFile = validParameter.validFile(parameters, "fasta", true);
			if (fastaFile == "not found") { cout << "fasta is a required parameter for the screen.seqs command." << endl; abort = true; }
			else if (fastaFile == "not open") { abort = true; }	
		
		
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
			else if(temp == "not open")	{	abort = 0;		}
			else						{	qFileName = temp;	}
			
			temp = validParameter.validFile(parameters, "qthreshold", false);	if (temp == "not found") { temp = "0"; }
			convert(temp, qThreshold);

			temp = validParameter.validFile(parameters, "qaverage", false);		if (temp == "not found") { temp = "0"; }
			convert(temp, qAverage);
			
			temp = validParameter.validFile(parameters, "allfiles", false);		if (temp == "not found") { temp = "F"; }
			allFiles = isTrue(temp);

			if(allFiles && oligoFile == ""){
				cout << "You selected allfiles, but didn't enter an oligos file.  Ignoring the allfiles request." << endl;
			}
			
			if(!flip && oligoFile=="" && !maxLength && !minLength && (maxAmbig==-1) && !maxHomoP && qFileName == ""){		
				cout << "You didn't set any options... quiting command." << endl;
				abort = true;
			}
			
		}

	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the TrimSeqsCommand class Function TrimSeqsCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the TrimSeqsCommand class function TrimSeqsCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}
//**********************************************************************************************************************

void TrimSeqsCommand::help(){
	try {
		cout << "The trim.seqs command reads a fastaFile and creates ....." << "\n";
		cout << "The trim.seqs command parameters are fasta, flip, oligos, maxambig, maxhomop, minlength and maxlength." << "\n";
		cout << "The fasta parameter is required." << "\n";
		cout << "The flip parameter .... The default is 0." << "\n";
		cout << "The oligos parameter .... The default is ""." << "\n";
		cout << "The maxambig parameter .... The default is -1." << "\n";
		cout << "The maxhomop parameter .... The default is 0." << "\n";
		cout << "The minlength parameter .... The default is 0." << "\n";
		cout << "The maxlength parameter .... The default is 0." << "\n";
		cout << "The trim.seqs command should be in the following format: " << "\n";
		cout << "trim.seqs(fasta=yourFastaFile, flip=yourFlip, oligos=yourOligos, maxambig=yourMaxambig,  " << "\n";
		cout << "maxhomop=yourMaxhomop, minlength=youMinlength, maxlength=yourMaxlength)  " << "\n";	
		cout << "Example trim.seqs(fasta=abrecovery.fasta, flip=..., oligos=..., maxambig=..., maxhomop=..., minlength=..., maxlength=...)." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFasta)." << "\n" << "\n";

	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the TrimSeqsCommand class Function help. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the TrimSeqsCommand class function help. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}


//***************************************************************************************************************

TrimSeqsCommand::~TrimSeqsCommand(){	/*	do nothing	*/	}

//***************************************************************************************************************

int TrimSeqsCommand::execute(){
	try{
	
		if (abort == true) { return 0; }
	
		vector<ofstream*> groupFileNames;
		vector<ofstream*> fastaFileNames;
		if(oligoFile != "")	{	getOligos(fastaFileNames, groupFileNames);	}

		ifstream inFASTA;
		openInputFile(fastaFile, inFASTA);

		ofstream outFASTA;
		string trimSeqFile = getRootName(fastaFile) + "trim.fasta";
		openOutputFile(trimSeqFile, outFASTA);
		
		ofstream outGroups;
		string groupFile = getRootName(fastaFile) + "groups"; 
		openOutputFile(groupFile, outGroups);
		
		ofstream scrapFASTA;
		string scrapSeqFile = getRootName(fastaFile) + "scrap.fasta";
		openOutputFile(scrapSeqFile, scrapFASTA);
		
		ifstream qFile;
		if(qFileName != "")	{	openInputFile(qFileName, qFile);	}
		
		bool success;
		
		while(!inFASTA.eof()){
			Sequence currSeq(inFASTA);
			string origSeq = currSeq.getUnaligned();
			int group;
			string trashCode = "";
			
			if(qFileName != ""){
				if(qThreshold != 0)		{	success = stripQualThreshold(currSeq, qFile);	}
				else if(qAverage != 0)	{	success = cullQualAverage(currSeq, qFile);		}
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
				if(!success){	trashCode += 'l';	}
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
				currSeq.printSequence(outFASTA);
				if(barcodes.size() != 0){
					outGroups << currSeq.getName() << '\t' << groupVector[group] << endl;
					
					if(allFiles){
						*groupFileNames[group] << currSeq.getName() << '\t' << groupVector[group] << endl;					
						currSeq.printSequence(*fastaFileNames[group]);					
					}
				}
			}
			else{
				currSeq.setName(currSeq.getName() + '|' + trashCode);
				currSeq.setUnaligned(origSeq);
				currSeq.printSequence(scrapFASTA);
			}
			gobble(inFASTA);
		}
		inFASTA.close();
		outFASTA.close();
		scrapFASTA.close();
		outGroups.close();
		
		for(int i=0;i<groupFileNames.size();i++){
			groupFileNames[i]->close();
			delete groupFileNames[i];

			fastaFileNames[i]->close();
			delete fastaFileNames[i];
		}
		
		return 0;		
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the TrimSeqsCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the TrimSeqsCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

//***************************************************************************************************************

void TrimSeqsCommand::getOligos(vector<ofstream*>& outFASTAVec, vector<ofstream*>& outGroupsVec){
	
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
				revPrimer.push_back(oligo);
			}
			else if(type == "barcode"){
				inOligos >> group;
				barcodes[oligo]=index++;
				groupVector.push_back(group);
					
				if(allFiles){
					outFASTAVec.push_back(new ofstream((getRootName(fastaFile) + group + ".fasta").c_str(), ios::ate));
					outGroupsVec.push_back(new ofstream((getRootName(fastaFile) + group + ".groups").c_str(), ios::ate));
				}
			}
		}
	}
	
	numFPrimers = forPrimer.size();
	numRPrimers = revPrimer.size();
}

//***************************************************************************************************************

bool TrimSeqsCommand::stripBarcode(Sequence& seq, int& group){
	
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

//***************************************************************************************************************

bool TrimSeqsCommand::stripForward(Sequence& seq){
	
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

//***************************************************************************************************************

bool TrimSeqsCommand::stripReverse(Sequence& seq){
	
	string rawSequence = seq.getUnaligned();
	bool success = 0;	//guilty until proven innocent
	
	for(int i=0;i<numRPrimers;i++){
		string oligo = revPrimer[i];
		
		if(rawSequence.length() < oligo.length()){
			success = 0;
			break;
		}
		
		if(compareDNASeq(oligo, rawSequence.substr(rawSequence.length()-oligo.length(),oligo.length()))){
			seq.setUnaligned(rawSequence.substr(rawSequence.length()-oligo.length()));
			success = 1;
			break;
		}
	}	
	return success;
	
}

//***************************************************************************************************************

bool TrimSeqsCommand::cullLength(Sequence& seq){
	
	int length = seq.getNumBases();
	bool success = 0;	//guilty until proven innocent
	
	if(length >= minLength && maxLength == 0)			{	success = 1;	}
	else if(length >= minLength && length <= maxLength)	{	success = 1;	}
	else												{	success = 0;	}
	
	return success;
	
}

//***************************************************************************************************************

bool TrimSeqsCommand::cullHomoP(Sequence& seq){
	
	int longHomoP = seq.getLongHomoPolymer();
	bool success = 0;	//guilty until proven innocent
	
	if(longHomoP <= maxHomoP){	success = 1;	}
	else					{	success = 0;	}
	
	return success;
	
}

//***************************************************************************************************************

bool TrimSeqsCommand::cullAmbigs(Sequence& seq){
	
	int numNs = seq.getAmbigBases();
	bool success = 0;	//guilty until proven innocent
	
	if(numNs <= maxAmbig)	{	success = 1;	}
	else					{	success = 0;	}
	
	return success;
	
}

//***************************************************************************************************************

bool TrimSeqsCommand::compareDNASeq(string oligo, string seq){
	
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

//***************************************************************************************************************

bool TrimSeqsCommand::stripQualThreshold(Sequence& seq, ifstream& qFile){
	
	string rawSequence = seq.getUnaligned();
	int seqLength = rawSequence.length();
	string name;
	
	qFile >> name;
	if(name.substr(1) != seq.getName())	{	cout << "sequence name mismatch btwn fasta and qual file" << endl;	}
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

//***************************************************************************************************************

bool TrimSeqsCommand::cullQualAverage(Sequence& seq, ifstream& qFile){
	
	string rawSequence = seq.getUnaligned();
	int seqLength = seq.getNumBases();
	bool success = 0;	//guilty until proven innocent
	string name;
	
	qFile >> name;
	if(name.substr(1) != seq.getName())	{	cout << "sequence name mismatch btwn fasta and qual file" << endl;	}
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

//***************************************************************************************************************


