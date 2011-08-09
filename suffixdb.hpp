#ifndef SUFFIXDB_HPP
#define SUFFIXDB_HPP

/*
 *  suffixdb.hpp
 *  
 *
 *  Created by Pat Schloss on 12/16/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 *	This is a child class of the Database abstract datatype.  The class is basically a database of suffix trees and an
 *	encapsulation of the method for finding the most similar tree to an inputted sequence.  the suffixForest object
 *	is a vector of SuffixTrees, with each template sequence being represented by a different SuffixTree.  The class also
 *	provides a method to take an unaligned sequence and find the closest sequence in the suffixForest.  The search
 *	method is inspired by the article and Perl source code provided at http://www.ddj.com/web-development/184416093.  I 
 *	would estimate that the time complexity is O(LN) for each search, which is slower than the kmer searching, but 
 *	faster than blast
 *
 */

#include "mothur.h"
#include "database.hpp"
#include "suffixtree.hpp"
//class SuffixTree;

class SuffixDB : public Database {
	
public:
	SuffixDB(int);
	SuffixDB();
	SuffixDB(const SuffixDB& sdb) : count(sdb.count), Database(sdb) {
		for (int i = 0; i < sdb.suffixForest.size(); i++) {
			SuffixTree temp(sdb.suffixForest[i]);
			suffixForest.push_back(temp);
		}
	}
	~SuffixDB();
	
	void generateDB() {}; //adding sequences generates the db
	void addSequence(Sequence);
	vector<int> findClosestSequences(Sequence*, int);

private:
	vector<SuffixTree> suffixForest;
	int count;
};

#endif
