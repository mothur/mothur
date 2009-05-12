#ifndef READNEXUS_H
#define READNEXUS_H

/*
 *  readnuxus.h
 *  Mothur
 *
 *  Created by Thomas Ryabin on 4/22/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */
using namespace std;

#include "readseqs.h"
#include "globaldata.hpp"
#include "sequencedb.h"
#include "mothur.h"

/**********************************************************************************/

class ReadNexus : public ReadSeqs {

	public:
		ReadNexus(string);
		~ReadNexus();
		void read();
		SequenceDB* getDB();				
};

#endif