/*
 *  suffixdb.cpp
 *  
 *
 *  Created by Pat Schloss on 12/16/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 *	This is a child class of the Database abstract datatype.  The class is basically a database of suffix trees and an
 *	encapsulation of the method for finding the most similar tree to an inputted sequence.  the suffixForest objecct
 *	is a vector of SuffixTrees, with each template sequence being represented by a different SuffixTree.  The class also
 *	provides a method to take an unaligned sequence and find the closest sequence in the suffixForest.  The search
 *	method is inspired by the article and Perl source code provided at http://www.ddj.com/web-development/184416093.  I 
 *	would estimate that the time complexity is O(LN) for each search, which is slower than the kmer searching, but 
 *	faster than blast
 *
 */

#include "database.hpp"
#include "sequence.hpp"
#include "suffixtree.hpp"
#include "suffixdb.hpp"

/**************************************************************************************************/

SuffixDB::SuffixDB(int numSeqs) : Database() {
	suffixForest.resize(numSeqs);
	count = 0;
}
/**************************************************************************************************/
//assumes sequences have been added using addSequence
vector<int> SuffixDB::findClosestSequences(Sequence* candidateSeq, int num){
	try {
		vector<int> topMatches;
		string processedSeq = candidateSeq->convert2ints();		//	the candidate sequence needs to be a string of ints
		
		vector<seqMatch> seqMatches;
		for(int i=0;i<suffixForest.size();i++){					//	scan through the forest and see what the minimum
			int count = suffixForest[i].countSuffixes(processedSeq);	//	return score is and keep track of the
			seqMatch temp(i, count);
			seqMatches.push_back(temp);
		}
		
		//sorts putting smallest matches first
		sort(seqMatches.begin(), seqMatches.end(), compareSeqMatchesReverse);
		
		searchScore = seqMatches[0].match;
		searchScore = 100 * (1. - searchScore / (float)processedSeq.length());
		
		//save top matches
		for (int i = 0; i < num; i++) {
			topMatches.push_back(seqMatches[i].seq);
		}

		//	return the Sequence object that has the minimum score
		return topMatches;
	}
	catch(exception& e) {
		errorOut(e, "SuffixDB", "findClosestSequences");
		exit(1);
	}	
}
/**************************************************************************************************/
//adding the sequences generates the db
void SuffixDB::addSequence(Sequence seq) {
	try {
		suffixForest[count].loadSequence(seq);		
		count++;
	}
	catch(exception& e) {
		errorOut(e, "SuffixDB", "addSequence");
		exit(1);
	}
}
/**************************************************************************************************/

SuffixDB::~SuffixDB(){														
	for (int i = (suffixForest.size()-1); i >= 0; i--) {  suffixForest.pop_back();  }
}
/**************************************************************************************************/
