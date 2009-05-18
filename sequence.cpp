/*
 *  sequence.cpp
 *  
 *
 *  Created by Pat Schloss on 12/15/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 */

using namespace std;

#include "sequence.hpp"

/***********************************************************************/

Sequence::Sequence()  {}

/***********************************************************************/

Sequence::Sequence(string newName, string sequence) {
	name = newName;
	if(sequence.find_first_of('-') != string::npos) {
		setAligned(sequence);
	}
	setUnaligned(sequence);
}
//********************************************************************************************************************

Sequence::Sequence(ifstream& fastaFile){
	
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

string Sequence::convert2ints() {
	
	if(unaligned == "")	{	/* need to throw an error */	}
	
	string processed;
	
	for(int i=0;i<unaligned.length();i++) {
		if(toupper(unaligned[i]) == 'A')			{	processed += '0';	}
		else if(toupper(unaligned[i]) == 'C')	{	processed += '1';	}
		else if(toupper(unaligned[i]) == 'G')	{	processed += '2';	}
		else if(toupper(unaligned[i]) == 'T')	{	processed += '3';	}
		else if(toupper(unaligned[i]) == 'U')	{	processed += '3';	}
		else									{	processed += '4';	}
	}
	return processed;
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
	
}

//********************************************************************************************************************

void Sequence::setAligned(string sequence){
	aligned = sequence;
}

//********************************************************************************************************************

void Sequence::setPairwise(string sequence){
	pairwise = sequence;
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

int Sequence::getLength(){
	if(unaligned.length() > aligned.length())
		return unaligned.length();
	return aligned.length();
}

//********************************************************************************************************************

void Sequence::printSequence(ostream& out){
	string toPrint = unaligned;
	if(aligned.length() > unaligned.length())
		toPrint = aligned;
	out << ">" << name << "\n" << toPrint << "\n";
}

//********************************************************************************************************************

int Sequence::getUnalignLength(){
	return unaligned.length();
}

//********************************************************************************************************************

int Sequence::getAlignLength(){
	return aligned.length();
}

//********************************************************************************************************************



