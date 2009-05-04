/*
 *  readphylip.cpp
 *  Mothur
 *
 *  Created by Thomas Ryabin on 4/24/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "readseqsphylip.h"
#include <iostream>
#include <fstream>

/*******************************************************************************/
bool ReadPhylip::isSeq(string seq) {
	string validChars[] = {"A","G","C","T","U","N","-"};
	
	for(int i = 0; i < seq.length(); i++) {
		bool valid = false;
		string c = seq.substr(i,1);
		for(int k = 0; k < 7; k++)
			if(c.compare(validChars[k]) == 0) {
				valid = true;
				k = 7;
			}
		if(!valid)
			return false;
	}
	
	return true;
}

/*******************************************************************************/
ReadPhylip::ReadPhylip(string file) {
	try {
		openInputFile(file, filehandle);
		phylipFile = file;
		globaldata = GlobalData::getInstance();
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ReadTree class Function ReadTree. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ReadTree class function ReadTree. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		
}
/*******************************************************************************/
ReadPhylip::~ReadPhylip(){
//	for(int i = 0; i < sequencedb.getNumSeqs(); i++)
//		delete sequencedb.get(i);
}
/*******************************************************************************/
void ReadPhylip::read() {
	string temp;
	string name;
	string sequence;
	
	int count = 0;
	int letterCount = 0;
	int numCols = 0;
	filehandle >> temp;
	int numSeqs = atoi(temp.c_str());
	filehandle >> temp;
	int numLetters = atoi(temp.c_str());
	
	bool firstDone = false;	
	bool last = false;
	filehandle >> name;
	
	while(!filehandle.eof()) {
		if(!firstDone) {
			sequence = "";
			if(count == 0) {
				filehandle >> temp;
				while(isSeq(temp)) {
					sequence += temp;
					numCols++;
					filehandle >> temp;
				}
				letterCount += sequence.length();
			}
			else {
				for(int i = 0; i < numCols; i++) {
					filehandle >> temp;
					sequence += temp;
				}
				if(count < numSeqs-1)
					filehandle >> temp;
			}
			Sequence newSeq(name, sequence);
			sequencedb.add(newSeq);
			if(count < numSeqs-1)
				name = temp;
		} 	
		else {
			sequence = "";
			for(int i = 0; i < numCols; i++) {
				filehandle >> temp;
				sequence += temp;
				if(count == 0)
					letterCount += temp.length();
				if(letterCount == numLetters && count == 0) {
					numCols = i + 1;
					i = numCols;
				}
			}
			if(!(last && count == 0))
				sequencedb.set(count, sequencedb.get(count).getAligned() + sequence);
			if(letterCount == numLetters && count == 0)
				last = true;
		}
		
		count++;
		
		if(count == numSeqs) {
			firstDone = true;
			count = 0;
		}
	}
	filehandle.close();
}

/*********************************************************************************/
SequenceDB* ReadPhylip::getDB() {
	return &sequencedb;
}
