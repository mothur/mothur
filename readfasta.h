#ifndef READFASTA_H
#define READFASTA_H

/*
 *  readfasta.h
 *  Mothur
 *
 *  Created by Thomas Ryabin on 4/21/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

using namespace std;

#include "readseqs.h"
#include "globaldata.hpp"
#include "sequencedb.h"
#include "mothur.h"

/**********************************************************************************/

class ReadFasta : public ReadSeqs {

	public:
		ReadFasta(string);
		~ReadFasta();
		void read();
		SequenceDB* getDB();		
};

#endif