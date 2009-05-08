#ifndef READPHYLIP_H
#define READPHYLIP_H

/*
 *  readphylip.h
 *  Mothur
 *
 *  Created by Thomas Ryabin on 4/24/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */
using namespace std;

#include "globaldata.hpp"
#include "sequencedb.h"
#include "mothur.h"

/**********************************************************************************/

class ReadPhylip {

	public:
		ReadPhylip(string);
		~ReadPhylip();
		void read();
		SequenceDB* getDB();		
	
	private:
		GlobalData* globaldata;
		string phylipFile;
		ifstream filehandle;
		SequenceDB sequencedb;
		int readOk; // readOk = 0 means success, readOk = 1 means error(s).
		
		bool isSeq(string);
};

#endif