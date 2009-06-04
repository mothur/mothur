/*
 *  readseqs.cpp
 *  Mothur
 *
 *  Created by Thomas Ryabin on 5/11/09.
 *  Copyright 2009Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "readseqs.h"

/*******************************************************************************/
ReadSeqs::ReadSeqs(string file) {
	try {
		openInputFile(file, filehandle);
		seqFile = file;
		globaldata = GlobalData::getInstance();
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ReadSeqs class Function ReadSeqs. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ReadSeqs class function ReadSeqs. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		
}
/*******************************************************************************/
ReadSeqs::~ReadSeqs(){
	//for(int i = 0; i < sequencedb.getNumSeqs(); i++)
		//delete sequencedb.get(i);
}
/*******************************************************************************/
void ReadSeqs::read() {
}

/*********************************************************************************/
SequenceDB* ReadSeqs::getDB() {
	return &sequencedb;
}
