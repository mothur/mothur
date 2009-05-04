#ifndef READNEXUS_H
#define READNEXUS_H

/*
 *  readnuxus.h
 *  Mothur
 *
 *  Created by Thomas Ryabin on 4/22/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */
using namespace std;

#include "globaldata.hpp"
#include "sequencedb.h"
#include "mothur.h"

/**********************************************************************************/

class ReadNexus {

	public:
		ReadNexus(string);
		~ReadNexus();
		void read();
		SequenceDB* getDB();		
	
	private:
		GlobalData* globaldata;
		string nexusFile;
		ifstream filehandle;
		SequenceDB sequencedb;
		int readOk; // readOk = 0 means success, readOk = 1 means error(s).
		
			
};

#endif