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
		
		globaldata = GlobalData::getInstance();
		abort = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string AlignArray[] =  {"fasta", "flip", "oligos", "maxambig", "maxhomop", "minlength", "maxlength"};
			vector<string> myArray (AlignArray, AlignArray+(sizeof(AlignArray)/sizeof(string)));
			
			parser = new OptionParser();
			parser->parse(option, parameters);   	delete parser; 
			
			ValidParameters* validParameter = new ValidParameters();
		
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter->isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//check for required parameters
			fastafile = validParameter->validFile(parameters, "fasta", true);
			if (fastafile == "not found") { cout << "fasta is a required parameter for the screen.seqs command." << endl; abort = true; }
			else if (fastafile == "not open") { abort = true; }	
			else { globaldata->setFastaFile(fastafile); }
		
		
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			string temp;
			temp = validParameter->validFile(parameters, "flip", false);			if (temp == "not found") { temp = "0"; }
			if(isTrue(temp))	{	flip = 1;	}
		
			temp = validParameter->validFile(parameters, "oligos", false);			if (temp == "not found") { temp = ""; }
			if(temp != "")		{	oligos = 1;	 } 
			else {  oligos = 0;	 }

			temp = validParameter->validFile(parameters, "maxambig", false);		if (temp == "not found") { temp = "-1"; }
			convert(temp, maxAmbig);  

			temp = validParameter->validFile(parameters, "maxhomop", false);		if (temp == "not found") { temp = "0"; }
			convert(temp, maxHomoP);  

			temp = validParameter->validFile(parameters, "minlength", false);		if (temp == "not found") { temp = "0"; }
			convert(temp, minLength); 
			
			temp = validParameter->validFile(parameters, "maxlength", false);		if (temp == "not found") { temp = "0"; }
			convert(temp, maxLength); 
			
			if(!flip && !oligos && !maxLength && !minLength && (maxAmbig==-1) && !maxHomoP ){	cout << "huh?" << endl;	}
			
			delete validParameter;
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
		cout << "The trim.seqs command reads a fastafile and creates ....." << "\n";
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
	
		getOligos();
		
		ifstream inFASTA;
		openInputFile(fastafile, inFASTA);

		ofstream outFASTA;
		string trimSeqFile = getRootName(fastafile) + "trim.fasta";
		openOutputFile(trimSeqFile, outFASTA);
		
		ofstream outGroups;
		string groupFile = getRootName(fastafile) + "groups"; 
		openOutputFile(groupFile, outGroups);

		ofstream scrapFASTA;
		string scrapSeqFile = getRootName(fastafile) + "scrap.fasta";
		openOutputFile(scrapSeqFile, scrapFASTA);

		bool success;
		
		while(!inFASTA.eof()){
			Sequence currSeq(inFASTA);
			string origSeq = currSeq.getUnaligned();
			string group;
			string trashCode = "";

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
				outGroups << currSeq.getName() << '\t' << group << endl;
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

void TrimSeqsCommand::getOligos(){

	ifstream inOligos;
	//openInputFile(globaldata->getOligosFile(), inOligos);

	string type, oligo, group;
	
	while(!inOligos.eof()){
		inOligos >> type;
		
		if(type == "forward"){
			inOligos >> oligo;
			forPrimer.push_back(oligo);
		}
		else if(type == "reverse"){
			inOligos >> oligo;
			revPrimer.push_back(oligo);
		}
		else if(type == "barcode"){
			inOligos >> oligo >> group;
			barcodes[oligo]=group;
		}
		else if(type[0] == '#'){
			char c;
			while ((c = inOligos.get()) != EOF)	{	if (c == 10){	break;	}	} // get rest of line
		}
		
		gobble(inOligos);
	}

	numFPrimers = forPrimer.size();
	numRPrimers = revPrimer.size();
}

//***************************************************************************************************************

bool TrimSeqsCommand::stripBarcode(Sequence& seq, string& group){
	
	string rawSequence = seq.getUnaligned();
	bool success = 0;	//guilty until proven innocent

	for(map<string,string>::iterator it=barcodes.begin();it!=barcodes.end();it++){
		string oligo = it->first;
		
		if(rawSequence.length() < oligo.length()){	//let's just assume that the barcodes are the same length
			success = 0;
			break;
		}
		
		if (rawSequence.compare(0,oligo.length(),oligo) == 0){
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
		
		if (rawSequence.compare(0,oligo.length(),oligo) == 0){
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
		
		if(rawSequence.compare(rawSequence.length()-oligo.length(),oligo.length(),oligo) == 0){
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
	
	if(numNs <= maxAmbig){	success = 1;	}
	else					{	success = 0;	}
	
	return success;
	
}

//***************************************************************************************************************
