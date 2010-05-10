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
	m = MothurOut::getInstance();
	initialize();
}
/***********************************************************************/
Sequence::Sequence(string newName, string sequence) {
	try {
		m = MothurOut::getInstance();
		initialize();	
		name = newName;
		
		//setUnaligned removes any gap characters for us
		setUnaligned(sequence);
		setAligned(sequence);
	}
	catch(exception& e) {
		m->errorOut(e, "Sequence", "Sequence");
		exit(1);
	}			
}
/***********************************************************************/
Sequence::Sequence(string newName, string sequence, string justUnAligned) {
	try {
		m = MothurOut::getInstance();
		initialize();	
		name = newName;
		
		//setUnaligned removes any gap characters for us
		setUnaligned(sequence);
	}
	catch(exception& e) {
		m->errorOut(e, "Sequence", "Sequence");
		exit(1);
	}			
}

//********************************************************************************************************************
//this function will jump over commented out sequences, but if the last sequence in a file is commented out it makes a blank seq
Sequence::Sequence(istringstream& fastaString){
	try {
		m = MothurOut::getInstance();
	
		initialize();
		fastaString >> name;
		name = name.substr(1);
		string sequence;
	
		//read comments
		while ((name[0] == '#') && fastaString) { 
			while (!fastaString.eof())	{	char c = fastaString.get(); if (c == 10 || c == 13){	break;	}	} // get rest of line if there's any crap there
			sequence = getCommentString(fastaString);
			
			if (fastaString) {  
				fastaString >> name;  
				name = name.substr(1);	
			}else { 
				name = "";
				break;
			}
		}
		
		while (!fastaString.eof())	{	char c = fastaString.get();  if (c == 10 || c == 13){	break;	}	} // get rest of line if there's any crap there
		
		sequence = getSequenceString(fastaString);		
		setAligned(sequence);	
		//setUnaligned removes any gap characters for us						
		setUnaligned(sequence);		
	}
	catch(exception& e) {
		m->errorOut(e, "Sequence", "Sequence");
		exit(1);
	}								
}
//********************************************************************************************************************
//this function will jump over commented out sequences, but if the last sequence in a file is commented out it makes a blank seq
Sequence::Sequence(istringstream& fastaString, string JustUnaligned){
	try {
		m = MothurOut::getInstance();
	
		initialize();
		fastaString >> name;
		name = name.substr(1);
		string sequence;
	
		//read comments
		while ((name[0] == '#') && fastaString) { 
			while (!fastaString.eof())	{	char c = fastaString.get(); if (c == 10 || c == 13){	break;	}	} // get rest of line if there's any crap there
			sequence = getCommentString(fastaString);
			
			if (fastaString) {  
				fastaString >> name;  
				name = name.substr(1);	
			}else { 
				name = "";
				break;
			}
		}
		
		while (!fastaString.eof())	{	char c = fastaString.get();  if (c == 10 || c == 13){	break;	}	} // get rest of line if there's any crap there
		
		sequence = getSequenceString(fastaString);		
		
		//setUnaligned removes any gap characters for us						
		setUnaligned(sequence);		
	}
	catch(exception& e) {
		m->errorOut(e, "Sequence", "Sequence");
		exit(1);
	}								
}


//********************************************************************************************************************
//this function will jump over commented out sequences, but if the last sequence in a file is commented out it makes a blank seq
Sequence::Sequence(ifstream& fastaFile){
	try {
		m = MothurOut::getInstance();
		initialize();
		fastaFile >> name;
		name = name.substr(1);
		string sequence;
		
		//read comments
		while ((name[0] == '#') && fastaFile) { 
			while (!fastaFile.eof())	{	char c = fastaFile.get(); if (c == 10 || c == 13){	break;	}	} // get rest of line if there's any crap there
			sequence = getCommentString(fastaFile);
			
			if (fastaFile) {  
				fastaFile >> name;  
				name = name.substr(1);	
			}else { 
				name = "";
				break;
			}
		}
		
		//read real sequence
		while (!fastaFile.eof())	{	char c = fastaFile.get(); if (c == 10 || c == 13){	break;	}	} // get rest of line if there's any crap there
		
		sequence = getSequenceString(fastaFile);		
		
		setAligned(sequence);	
		//setUnaligned removes any gap characters for us						
		setUnaligned(sequence);	
	}
	catch(exception& e) {
		m->errorOut(e, "Sequence", "Sequence");
		exit(1);
	}							
}
//********************************************************************************************************************
//this function will jump over commented out sequences, but if the last sequence in a file is commented out it makes a blank seq
Sequence::Sequence(ifstream& fastaFile, string JustUnaligned){
	try {
		m = MothurOut::getInstance();
		initialize();
		fastaFile >> name;
		name = name.substr(1);
		string sequence;
		
		//read comments
		while ((name[0] == '#') && fastaFile) { 
			while (!fastaFile.eof())	{	char c = fastaFile.get(); if (c == 10 || c == 13){	break;	}	} // get rest of line if there's any crap there
			sequence = getCommentString(fastaFile);
			
			if (fastaFile) {  
				fastaFile >> name;  
				name = name.substr(1);	
			}else { 
				name = "";
				break;
			}
		}
		
		//read real sequence
		while (!fastaFile.eof())	{	char c = fastaFile.get(); if (c == 10 || c == 13){	break;	}	} // get rest of line if there's any crap there
		
		sequence = getSequenceString(fastaFile);		
		
		//setUnaligned removes any gap characters for us						
		setUnaligned(sequence);	
	}
	catch(exception& e) {
		m->errorOut(e, "Sequence", "Sequence");
		exit(1);
	}							
}

//********************************************************************************************************************
string Sequence::getSequenceString(ifstream& fastaFile) {
	try {
		char letter;
		string sequence = "";	
		
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
		
		return sequence;
	}
	catch(exception& e) {
		m->errorOut(e, "Sequence", "getSequenceString");
		exit(1);
	}
}
//********************************************************************************************************************
//comment can contain '>' so we need to account for that
string Sequence::getCommentString(ifstream& fastaFile) {
	try {
		char letter;
		string sequence = "";
		
		while(fastaFile){
			letter=fastaFile.get();
			if((letter == '\r') || (letter == '\n')){  
				gobble(fastaFile);  //in case its a \r\n situation
				break;
			}
		}
		
		return sequence;
	}
	catch(exception& e) {
		m->errorOut(e, "Sequence", "getCommentString");
		exit(1);
	}
}
//********************************************************************************************************************
string Sequence::getSequenceString(istringstream& fastaFile) {
	try {
		char letter;
		string sequence = "";	
		
		while(!fastaFile.eof()){
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
		
		return sequence;
	}
	catch(exception& e) {
		m->errorOut(e, "Sequence", "getSequenceString");
		exit(1);
	}
}
//********************************************************************************************************************
//comment can contain '>' so we need to account for that
string Sequence::getCommentString(istringstream& fastaFile) {
	try {
		char letter;
		string sequence = "";
		
		while(fastaFile){
			letter=fastaFile.get();
			if((letter == '\r') || (letter == '\n')){  
				gobble(fastaFile);  //in case its a \r\n situation
				break;
			}
		}
		
		return sequence;
	}
	catch(exception& e) {
		m->errorOut(e, "Sequence", "getCommentString");
		exit(1);
	}
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
	
	if(sequence.find_first_of('.') != string::npos || sequence.find_first_of('-') != string::npos) {
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
	setUnaligned(sequence);	

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
	if(isAligned == 0)	{ return unaligned; }
	else				{  return aligned;  }
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
	aligned = temp;
	
}
/**************************************************************************************************/
