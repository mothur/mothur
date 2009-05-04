/*
 *  readnexusaln.cpp
 *  Mothur
 *
 *  Created by Thomas Ryabin on 4/22/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "readnexus.h"
#include <iostream>
#include <fstream>

/*******************************************************************************/
ReadNexus::ReadNexus(string file) {
	try {
		openInputFile(file, filehandle);
		nexusFile = file;
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
ReadNexus::~ReadNexus(){
//	for(int i = 0; i < sequencedb.getNumSeqs(); i++)
//		delete sequencedb.get(i);
}
/*******************************************************************************/
void ReadNexus::read() {
	string temp;
	string name;
	string sequence;
	for(int i = 0; i < 5; i++)
		filehandle >> temp;
	
	int numSeqs = atoi(temp.substr(temp.find_first_of("=")+1, temp.length() - temp.find_first_of("=") - 1).c_str());
	
	for(int i = 0; i < 9; i++)
		filehandle >> temp;
	
	int count = 0;
	bool firstDone = false;	
	while(!filehandle.eof()){
		filehandle >> name;
		filehandle >> sequence;
		if(name.compare(";") != 0) {
			if(!firstDone) {
				Sequence newSeq(name, sequence);
				sequencedb.add(newSeq);
			} 
			else 
				sequencedb.set(count, sequencedb.get(count).getAligned() + sequence);
			
			count++;
			if(count == numSeqs) {
				if(!firstDone)
					firstDone = true;
				count = 0;
			}
		}
	}
	filehandle.close();
}

/*********************************************************************************/
SequenceDB* ReadNexus::getDB() {
	return &sequencedb;
}
