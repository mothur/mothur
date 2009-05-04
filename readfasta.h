#ifndef READFASTA_H
#define READFASTA_H

/*
 *  readfasta.h
 *  Mothur
 *
 *  Created by Thomas Ryabin on 4/21/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

using namespace std;

#include "globaldata.hpp"
#include "sequencedb.h"
#include "mothur.h"

/**********************************************************************************/

class ReadFasta {

	public:
		ReadFasta(string);
		~ReadFasta();
		void read();
		SequenceDB* getDB();		
	
	private:
		GlobalData* globaldata;
		string fastaFile;
		ifstream filehandle;
		SequenceDB sequencedb;
		int readOk; // readOk = 0 means success, readOk = 1 means error(s).
		
			
};

#endif