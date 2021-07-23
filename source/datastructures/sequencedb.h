#ifndef SEQUENCEDB_H
#define SEQUENCEDB_H

/*
 *  sequencedb.h
 *  Mothur
 *
 *  Created by Thomas Ryabin on 4/13/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */


/* This class is a container to store the sequences. */

#include "storagedatabase.hpp"
#include "sequence.hpp"

class SequenceDB : public StorageDatabase {
	
public:
	SequenceDB();
	SequenceDB(int);           //makes data that size
	SequenceDB(ifstream&);	   //reads file to fill data
    SequenceDB(ifstream&, int, vector< vector< int > >&, vector< int >&); //filehandle, kmersize, kmerdb, lengths
	SequenceDB(const SequenceDB& sdb) : data(sdb.data) {};
    SequenceDB(const SequenceDB& sdb, set<string> names); //creates a new sequenceDB containing only the reads in names 
	~SequenceDB();             //loops through data and delete each sequence

	int getNumSeqs();
    Sequence getSeq(int);         //returns sequence at that location
	void push_back(Sequence);  //adds unaligned sequence
	bool sameLength() { return samelength; }
		
private:
	vector<Sequence> data;
	
	

};

#endif
