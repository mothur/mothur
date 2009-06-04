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

SuffixDB::SuffixDB(string fastaFileName) : Database(fastaFileName) {

	suffixForest.resize(numSeqs);
	cout << "Generating the suffix tree database...\t";	cout.flush();
	for(int i=0;i<numSeqs;i++){								//	The parent class' constructor generates the vector of
		suffixForest[i].loadSequence(templateSequences[i]);	//	template Sequence objects.  Here each of these objects
	}														//	is used to generate a suffix tree, aka the suffix forest
	cout << "DONE." << endl << endl;	cout.flush();

}

/**************************************************************************************************/

Sequence* SuffixDB::findClosestSequence(Sequence* candidateSeq){

	int minValue = 2000;
	int closestSequenceNo = 0;
	string processedSeq = candidateSeq->convert2ints();		//	the candidate sequence needs to be a string of ints
	for(int i=0;i<suffixForest.size();i++){					//	scan through the forest and see what the minimum
		int count = suffixForest[i].countSuffixes(processedSeq, minValue);	//	return score is and keep track of the
		if(count == minValue){								//	template sequence index that corresponds to that score
			closestSequenceNo = i;
		}
	}
	searchScore = 100 * (1. - minValue / (float)processedSeq.length());
	return templateSequences[closestSequenceNo];			//	return the Sequence object that has the minimum score
	
}

/**************************************************************************************************/
