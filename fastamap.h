#ifndef FASTAMAP_H
#define FASTAMAP_H

/*
 *  fastamap.h
 *  mothur
 *
 *  Created by Sarah Westcott on 1/16/09.
 *  Copyright 2009 Schloss Lab UMASS AMherst. All rights reserved.
 *
 */
 
#include "mothur.h"


/* This class represents the fasta file.  It reads a fasta file a populates the internal data structure "data".
Data is a map where the key is the sequence and the value is a struct containing the sequences groupname, 
a list of the sequences names who have the same sequence and a number of how many sequence names there are. */


class FastaMap  {

public:
	FastaMap() {};
	~FastaMap() {};
	
	string getGroupName(string);  //pass a sequence name get its group
	int getGroupNumber(string);  //pass a sequence name get number of sequence in its group
	string getNames(string);	//pass a sequence get the string of names in the group separated by ','s.
	void push_back(string, string); //sequencename, sequence
	int sizeUnique();					//returns number of unique sequences
	void printNamesFile(ostream&);		//produces a 2 column file with the groupname in the first column and the names in the second column - a names file.
	void printCondensedFasta(ostream&);		//produces a fasta file.
	void readFastaFile(ifstream&);
	string getSequence(string);		//pass it a name of a sequence, it returns the sequence.

private:
	struct group {
		string groupname;					//the group name for identical sequences, will be set to the first sequence found.
		int groupnumber;					//the number of sequence names with the same sequence.
		string names;						//the names of the sequence separated by ','.
	};

	map<string, group>  data;  //sequence, groupinfo	- condensed representation of file
	map<string, string>  seqmap;  //name, sequence  -  uncondensed representation of file
	map<string, group>::iterator it;
	map<string, string>::iterator it2;
	
	string readName(ifstream&);
	string readSequence(ifstream&);
};

#endif
