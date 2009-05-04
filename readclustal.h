#ifndef READCLUSTAL_H
#define READCLUSTAL_H

/*
 *  readclustal.h
 *  Mothur
 *
 *  Created by Thomas Ryabin on 4/24/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */
using namespace std;

#include "globaldata.hpp"
#include "sequencedb.h"
#include "mothur.h"

/**********************************************************************************/

class ReadClustal {

	public:
		ReadClustal(string);
		~ReadClustal();
		void read();
		SequenceDB* getDB();		
	
	private:
		GlobalData* globaldata;
		string clustalFile;
		ifstream filehandle;
		SequenceDB sequencedb;
		int readOk; // readOk = 0 means success, readOk = 1 means error(s).
		
			
};

#endif