#ifndef READCLUSTAL_H
#define READCLUSTAL_H

/*
 *  readclustal.h
 *  Mothur
 *
 *  Created by Thomas Ryabin on 4/24/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "readseqs.h"
#include "globaldata.hpp"
#include "sequencedb.h"
#include "mothur.h"

/**********************************************************************************/

class ReadClustal : public ReadSeqs {

	public:
		ReadClustal(string);
		~ReadClustal();
		void read();
		SequenceDB* getDB();			
};

#endif