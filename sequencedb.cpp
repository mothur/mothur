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

SequenceDB::SequenceDB() {}

/***********************************************************************/

SequenceDB::SequenceDB(int newSize) {
	data.resize(newSize);
}

/***********************************************************************/

SequenceDB::SequenceDB(ifstream&) {}
	
/***********************************************************************/

int SequenceDB::getNumSeqs() {
	return data.size();
}

/***********************************************************************/

void SequenceDB::set(int index, string newUnaligned) {
	Sequence newSeq(data[index].getName(), newUnaligned);
	data[index] = newSeq;
}

/***********************************************************************/

void SequenceDB::set(int index, Sequence newSeq) {
	data[index] = newSeq;
}

/***********************************************************************/

Sequence SequenceDB::get(int index) {
	return data[index];
}

/***********************************************************************/

void SequenceDB::changeSize(int newSize) {
	data.resize(newSize);
}

/***********************************************************************/

void SequenceDB::clear() {
	data.clear();
}

/***********************************************************************/

int SequenceDB::size() {
	return data.size();
}

/***********************************************************************/

void SequenceDB::print(ofstream& out) {
	for(int i = 0; i < data.size(); i++)
		data[i].printSequence(out);
}
	
/***********************************************************************/

void SequenceDB::add(Sequence newSequence) {
	try {
		data.push_back(newSequence);
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the RAbundVector class Function push_back. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the RAbundVector class function push_back. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}