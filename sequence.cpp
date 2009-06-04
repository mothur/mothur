/*
 *  sequence.cpp
 *  
 *
 *  Created by Pat Schloss on 12/15/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 */

#include "sequence.hpp"

/***********************************************************************/

Sequence::Sequence(){
	initialize();
}

/***********************************************************************/

Sequence::Sequence(string newName, string sequence) {

	initialize();	
	name = newName;
	if(sequence.find_first_of('-') != string::npos) {
		setAligned(sequence);
	}
	setUnaligned(sequence);
	
}
//********************************************************************************************************************

Sequence::Sequence(ifstream& fastaFile){
	initialize();
	
	string accession;				//	provided a file handle to a fasta-formatted sequence file, read in the next
	fastaFile >> accession;			//	accession number and sequence we find...
	setName(accession);

	char letter;
	string sequence;
	
	while(fastaFile){
		letter= fastaFile.get();
		if(letter == '>'){
			fastaFile.putback(letter);
			break;
		}
		else if(isprint(letter)){
			letter = toupper(letter);
			if(letter == 'U'){letter = 'T';}
			sequence += letter;
		}
		
	}

	if(sequence.find_first_of('-') != string::npos){	//	if there are any gaps in the sequence, assume that it is
		setAligned(sequence);							//	an alignment file
	}
	setUnaligned(sequence);								//	also set the unaligned sequence file
}

//********************************************************************************************************************

void Sequence::initialize(){
	
	name = "";
	unaligned = "";
	aligned = "";
	pairwise = "";
	
	numBases = 0;
	alignmentLength = 0;
	isAligned = 0;
	startPos = -1;
	endPos = -1;
	longHomoPolymer = -1;
	ambigBases = -1;
	
}	

//********************************************************************************************************************

void Sequence::setName(string seqName) {
	if(seqName[0] == '>')	{	name = seqName.substr(1);	}
	else					{	name = seqName;				}
}

//********************************************************************************************************************

void Sequence::setUnaligned(string sequence){
	
	if(sequence.find_first_of('-') != string::npos) {
		string temp = "";
		for(int j=0;j<sequence.length();j++) {
			if(isalpha(sequence[j]))	{	temp += sequence[j];	}
		}
		unaligned = temp;
	}
	else {
		unaligned = sequence;
	}
	numBases = unaligned.length();
	
}

//********************************************************************************************************************

void Sequence::setAligned(string sequence){
	
	//if the alignment starts or ends with a gap, replace it with a period to indicate missing data
	aligned = sequence;
	alignmentLength = aligned.length();

	if(aligned[0] == '-'){
		for(int i=0;i<alignmentLength;i++){
			if(aligned[i] == '-'){
				aligned[i] = '.';
			}
			else{
				break;
			}
		}
		for(int i=alignmentLength-1;i>=0;i--){
			if(aligned[i] == '-'){
				aligned[i] = '.';
			}
			else{
				break;
			}
		}
	}
	isAligned = 1;	
}

//********************************************************************************************************************

void Sequence::setPairwise(string sequence){
	pairwise = sequence;
}

//********************************************************************************************************************

string Sequence::convert2ints() {
	
	if(unaligned == "")	{	/* need to throw an error */	}
	
	string processed;
	
	for(int i=0;i<unaligned.length();i++) {
		if(toupper(unaligned[i]) == 'A')		{	processed += '0';	}
		else if(toupper(unaligned[i]) == 'C')	{	processed += '1';	}
		else if(toupper(unaligned[i]) == 'G')	{	processed += '2';	}
		else if(toupper(unaligned[i]) == 'T')	{	processed += '3';	}
		else if(toupper(unaligned[i]) == 'U')	{	processed += '3';	}
		else									{	processed += '4';	}
	}
	return processed;
}

//********************************************************************************************************************

string Sequence::getName(){
	return name;
}

//********************************************************************************************************************

string Sequence::getAligned(){
	return aligned;
}

//********************************************************************************************************************

string Sequence::getPairwise(){
	return pairwise;
}

//********************************************************************************************************************

string Sequence::getUnaligned(){
	return unaligned;
}

//********************************************************************************************************************

int Sequence::getNumBases(){
	return numBases;
}

//********************************************************************************************************************

void Sequence::printSequence(ostream& out){

	out << ">" << name << endl;
	if(isAligned){
		out << aligned << endl;
	}
	else{
		out << unaligned << endl;
	}
}

//********************************************************************************************************************

int Sequence::getAlignLength(){
	return alignmentLength;
}

//********************************************************************************************************************

int Sequence::getAmbigBases(){
	if(ambigBases == -1){
		ambigBases = 0;
		for(int j=0;j<numBases;j++){
			if(unaligned[j] != 'A' && unaligned[j] != 'T' && unaligned[j] != 'G' && unaligned[j] != 'C'){
				ambigBases++;
			}
		}
	}	
	
	return ambigBases;
}

//********************************************************************************************************************

int Sequence::getLongHomoPolymer(){
	if(longHomoPolymer == -1){
		longHomoPolymer = 1;
		int homoPolymer = 1;
		for(int j=1;j<numBases;j++){
			if(unaligned[j] == unaligned[j-1]){
				homoPolymer++;
			}
			else{
				if(homoPolymer > longHomoPolymer){	longHomoPolymer = homoPolymer;	}
				homoPolymer = 1;
			}
		}
		if(homoPolymer > longHomoPolymer){	longHomoPolymer = homoPolymer;	}
	}
	return longHomoPolymer;
}

//********************************************************************************************************************

int Sequence::getStartPos(){
	if(endPos == -1){
		for(int j = 0; j < alignmentLength; j++) {
			if(aligned[j] != '.'){
				startPos = j + 1;
				break;
			}
		}
	}
	if(isAligned == 0){	startPos = 1;	}

	return startPos;
}

//********************************************************************************************************************

int Sequence::getEndPos(){
	if(endPos == -1){
		for(int j=alignmentLength-1;j>=0;j--){
			if(aligned[j] != '.'){
				endPos = j + 1;
				break;
			}
		}
	}
	if(isAligned == 0){	endPos = numBases;	}
	
	return endPos;
}

//********************************************************************************************************************

bool Sequence::getIsAligned(){
	return isAligned;
}

//********************************************************************************************************************

void Sequence::reverseComplement(){

	string temp;
	for(int i=numBases-1;i>=0;i--){
		if(unaligned[i] == 'A')		{	temp += 'T';	}
		else if(unaligned[i] == 'T'){	temp += 'A';	}
		else if(unaligned[i] == 'G'){	temp += 'C';	}
		else if(unaligned[i] == 'C'){	temp += 'G';	}
		else						{	temp += 'N';	}
	}
	unaligned = temp;
	
}

//********************************************************************************************************************
