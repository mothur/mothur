#ifndef READNEXUSALN_H
#define READNEXUSALN_H

/*
 *  readnexusaln.h
 *  Mothur
 *
 *  Created by Thomas Ryabin on 4/22/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

using namespace std;

#include "globaldata.hpp"
#include "sequencedb.h"
#include "utilities.hpp"

/**********************************************************************************/

class ReadNexus {

	public:
		ReadNexus(string);
		~ReadNexus();
		void read();		
	
	private:
		GlobalData* globaldata;
		string nexusFile;
		ifstream filehandle;
		SequenceDB sequencedb;
		int readOk; // readOk = 0 means success, readOk = 1 means error(s).
		
			
};

#endif