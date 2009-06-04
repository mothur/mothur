#ifndef READSEQS_H
#define READSEQS_H

/*
 *  readseqs.h
 *  Mothur
 *
 *  Created by Thomas Ryabin on 5/11/09.
 *  Copyright 2009Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "globaldata.hpp"
#include "sequencedb.h"
#include "mothur.h"

/**********************************************************************************/

class ReadSeqs {

	public:
		ReadSeqs(string);
		~ReadSeqs();
		virtual void read();
		virtual SequenceDB* getDB();		
	
	protected:
		GlobalData* globaldata;
		string seqFile;
		ifstream filehandle;
		SequenceDB sequencedb;
		int readOk; // readOk = 0 means success, readOk = 1 means error(s).		
};

#endif
