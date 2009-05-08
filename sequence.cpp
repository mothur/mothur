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

void Sequence::printSequence(ofstream& out){
	string toPrint = unaligned;
	if(aligned.length() > unaligned.length())
		toPrint = aligned;
	out << ">" << name << "\n" << toPrint << "\n";
}

//********************************************************************************************************************
