/*
 *  readclustal.cpp
 *  Mothur
 *
 *  Created by Thomas Ryabin on 4/24/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "readclustal.h"
#include <iostream>
#include <fstream>

/*******************************************************************************/
ReadClustal::ReadClustal(string file) {
	try {
		openInputFile(file, filehandle);
		clustalFile = file;
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
ReadClustal::~ReadClustal(){
//	for(int i = 0; i < sequencedb.getNumSeqs(); i++)
//		delete sequencedb.get(i);
}
/*******************************************************************************/
void ReadClustal::read() {
	string temp;
	string name;
	string sequence;
	string firstName = "";
	for(int i = 0; i < 6; i++)
		filehandle >> temp;
	
	int count = 0;
	int numSeqs = 0;
	bool firstDone = false;	
	
	while(!filehandle.eof()) {
		filehandle >> name;
		if(numSeqs != 0) {
			if(count == numSeqs)
				count = 0;
		}
		else if(!firstDone && firstName.compare("") == 0)
			firstName = name;
		else if(!firstDone && firstName.compare(name) == 0) {
			numSeqs = count;
			firstDone = true;
			count = 0;
		}

		if(name.find_first_of("*") == -1) {
			filehandle >> sequence;
			if(!firstDone) {
				Sequence newSeq(name, sequence);
				sequencedb.add(newSeq);
			} 
			else 
				sequencedb.set(count, sequencedb.get(count).getAligned() + sequence);
				
			count++;
		}
	}
	filehandle.close();
}

/*********************************************************************************/
SequenceDB* ReadClustal::getDB() {
	return &sequencedb;
}
