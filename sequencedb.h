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


using namespace std;

#include "sequence.hpp"
#include "calculator.h"




class SequenceDB {
	
public:
	SequenceDB();
	SequenceDB(int);           //makes data that size
	SequenceDB(ifstream&);	   //reads file to fill data
//	~SequenceDB();             //loops through data and delete each sequence

	int getNumSeqs();
	
	void set(int, string);     //unaligned - should also set length
	void set(int, Sequence);   //unaligned - should also set length
	Sequence get(int);         //returns sequence name at that location
	void add(Sequence);        //adds unaligned sequence
	void changeSize(int);      //resizes data
	void clear();              //clears data - remeber to loop through and delete the sequences inside or you will have a memory leak
	int size();                //returns datas size
	void print(ofstream&);      //loops through data using sequence class print
		
private:
	vector<Sequence> data;

};

#endif
