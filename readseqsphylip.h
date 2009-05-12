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

#include "readseqs.h"
#include "globaldata.hpp"
#include "sequencedb.h"
#include "mothur.h"

/**********************************************************************************/

class ReadPhylip : public ReadSeqs {

	public:
		ReadPhylip(string);
		~ReadPhylip();
		void read();
		SequenceDB* getDB();		
	
	private:
		bool isSeq(string);
};

#endif