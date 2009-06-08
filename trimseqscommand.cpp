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

TrimSeqsCommand::TrimSeqsCommand(){
	try {
		
		globaldata = GlobalData::getInstance();
		
		oligos = 0;		

		if(globaldata->getFastaFile() == ""){
			cout << "you need to at least enter a fasta file name" << endl;
		}
		
		if(isTrue(globaldata->getFlip()))			{	flip = 1;													}
		if(globaldata->getOligosFile() != "")		{	oligos = 1;													}
		
		if(globaldata->getMaxAmbig() != "-1")		{	maxAmbig = atoi(globaldata->getMaxAmbig().c_str());			}
		else										{	maxAmbig = -1;												}
		
		if(globaldata->getMaxHomoPolymer() != "-1")	{	maxHomoP = atoi(globaldata->getMaxHomoPolymer().c_str());	}
		else										{	maxHomoP = 0;												}
		
		if(globaldata->getMinLength() != "-1")		{	minLength = atoi(globaldata->getMinLength().c_str());		}
		else										{	minLength = 0;												}
		
		if(globaldata->getMaxLength() != "-1")		{	maxLength = atoi(globaldata->getMaxLength().c_str());		}
		else										{	maxLength = 0;												}
		
		if(!flip && !oligos && !maxLength && !minLength && (maxAmbig==-1) && !maxHomoP ){	cout << "huh?" << endl;	}

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

//***************************************************************************************************************

TrimSeqsCommand::~TrimSeqsCommand(){	/*	do nothing	*/	}

//***************************************************************************************************************

int TrimSeqsCommand::execute(){
	try{
		getOligos();
		
		ifstream inFASTA;
		openInputFile(globaldata->getFastaFile(), inFASTA);

		ofstream outFASTA;
		string trimSeqFile = getRootName(globaldata->getFastaFile()) + "trim.fasta";
		openOutputFile(trimSeqFile, outFASTA);
		
		ofstream outGroups;
		string groupFile = getRootName(globaldata->getFastaFile()) + "groups"; 
		openOutputFile(groupFile, outGroups);

		ofstream scrapFASTA;
		string scrapSeqFile = getRootName(globaldata->getFastaFile()) + "scrap.fasta";
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
	openInputFile(globaldata->getOligosFile(), inOligos);

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
