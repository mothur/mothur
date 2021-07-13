/*
 *  sequencedb.cpp
 *  Mothur
 *
 *  Created by Thomas Ryabin on 4/13/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "sequencedb.h"

/***********************************************************************/

SequenceDB::SequenceDB() : StorageDatabase() {}
/***********************************************************************/
//the clear function free's the memory
SequenceDB::~SequenceDB() { data.clear(); }

/***********************************************************************/

SequenceDB::SequenceDB(int newSize) : StorageDatabase() { data.resize(newSize, Sequence()); }

/***********************************************************************/

SequenceDB::SequenceDB(ifstream& filehandle) : StorageDatabase() {
	try{
       
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

SequenceDB::SequenceDB(const SequenceDB& sdb, set<string> names) : StorageDatabase() {
    try{
       
        int numSeqs = sdb.data.size();
        
        for (int i = 0; i < numSeqs; i++) {
            
            Sequence seqI = sdb.data[i];
            if (names.count(seqI.getName()) != 0) {
                if (length == 0) { length = seqI.getAligned().length(); }
                if (length != seqI.getAligned().length()) { samelength = false;  }
                data.push_back(seqI);
            }
        }
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
Sequence SequenceDB::getSeq(int index) {
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

