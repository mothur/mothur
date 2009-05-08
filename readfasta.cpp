/*
 *  readfasta.cpp
 *  Mothur
 *
 *  Created by Thomas Ryabin on 4/21/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "readfasta.h"
#include <iostream>
#include <fstream>

/*******************************************************************************/
ReadFasta::ReadFasta(string file) {
	try {
		openInputFile(file, filehandle);
		fastaFile = file;
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
ReadFasta::~ReadFasta(){
	//for(int i = 0; i < sequencedb.getNumSeqs(); i++)
		//delete sequencedb.get(i);
}
/*******************************************************************************/
void ReadFasta::read() {
	string name = "";
	string sequence = "";
	string temp;
	int count = 0;
	while(!filehandle.eof()){
		if(count == 0)
			filehandle >> temp;
		if(temp.substr(0,1).compare(">") == 0) {
			if(count != 0) {
				Sequence newSequence(name, sequence);
				sequencedb.add(newSequence);
				sequence = "";
			}
			else
				count++;
			name = temp.substr(1,temp.length()-1);
		}
		else 
			sequence += temp;
		
		filehandle >> temp;
	}
	Sequence newSequence(name, sequence);
	sequencedb.add(newSequence);

	filehandle.close();
}

/*********************************************************************************/
SequenceDB* ReadFasta::getDB() {
	return &sequencedb;
}
