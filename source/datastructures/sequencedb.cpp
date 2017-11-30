/*
 *  sequencedb.cpp
 *  Mothur
 *
 *  Created by Thomas Ryabin on 4/13/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "sequencedb.h"
#include "sequence.hpp"
#include "mothur.h"
#include "calculator.h"


/***********************************************************************/

SequenceDB::SequenceDB() {  m = MothurOut::getInstance();  length = 0; samelength = true; }
/***********************************************************************/
//the clear function free's the memory
SequenceDB::~SequenceDB() { data.clear(); }

/***********************************************************************/

SequenceDB::SequenceDB(int newSize) {
	data.resize(newSize, Sequence());
	length = 0; samelength = true;
}

/***********************************************************************/

SequenceDB::SequenceDB(ifstream& filehandle) {
	try{
		length = 0; samelength = true;
        Utils util;
		//read through file
		while (!filehandle.eof()) {
			//input sequence info into sequencedb
			Sequence newSequence(filehandle);
			
			if (newSequence.getName() != "") {   
				if (length == 0) { length = newSequence.getAligned().length(); }
                if (length != newSequence.getAligned().length()) { samelength = false;  }
				data.push_back(newSequence);  
			}
			
			//takes care of white space
			util.gobble(filehandle);
		}

		filehandle.close();
		
	}
	catch(exception& e) {
		m->errorOut(e, "SequenceDB", "SequenceDB");
		exit(1);
	}
}
/***********************************************************************/

int SequenceDB::getNumSeqs() {
	return data.size();
}

/***********************************************************************/
Sequence SequenceDB::get(int index) {
	return data[index];
}
/***********************************************************************/

void SequenceDB::push_back(Sequence newSequence) {
	try {
		if (length == 0) { length = newSequence.getAligned().length(); }
		if (length != newSequence.getAligned().length()) { samelength = false; }

		data.push_back(newSequence);
	}
	catch(exception& e) {
		m->errorOut(e, "SequenceDB", "push_back");
		exit(1);
	}
}

/***********************************************************************/

