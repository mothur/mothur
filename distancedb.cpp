/*
 *  distancedb.cpp
 *  
 *
 *  Created by Pat Schloss on 12/29/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 */


#include "database.hpp"
#include "sequence.hpp"
#include "distancedb.hpp"
#include "eachgapdist.h"

/**************************************************************************************************/
DistanceDB::DistanceDB() { 
	try {
		templateAligned = true;  
		templateSeqsLength = 0; 
		distCalculator = new eachGapDist();
	}
	catch(exception& e) {
		m->errorOut(e, "DistanceDB", "DistanceDB");
		exit(1);
	}	
}
/**************************************************************************************************/
void DistanceDB::addSequence(Sequence seq) {
	try {
		//are the template sequences aligned
		if (!isAligned(seq.getAligned())) { templateAligned = false; m->mothurOut(seq.getName() + " is not aligned. Sequences must be aligned to use the distance method."); m->mothurOutEndLine(); }
		
		if (templateSeqsLength == 0) { templateSeqsLength = seq.getAligned().length(); }
				
		data.push_back(seq);
	}
	catch(exception& e) {
		m->errorOut(e, "DistanceDB", "addSequence");
		exit(1);
	}	
}
/**************************************************************************************************/
//returns indexes to top matches
vector<int> DistanceDB::findClosestSequences(Sequence* query, int numWanted){
	try {
		vector<int> topMatches;
		bool templateSameLength = true;
		string sequence = query->getAligned();
		vector<seqDist> dists;
	
		if (numWanted > data.size()) { m->mothurOut("numwanted is larger than the number of template sequences, using "+ toString(data.size()) + "."); m->mothurOutEndLine(); numWanted = data.size(); }
		
		if (sequence.length() != templateSeqsLength) { templateSameLength = false; }
		
		if (templateSameLength && templateAligned) {
			//calc distance from this sequence to every sequence in the template
			for (int i = 0; i < data.size(); i++) {
				distCalculator->calcDist(*query, data[i]);
				float dist = distCalculator->getDist();
				
				//save distance to each template sequence
				seqDist temp(-1, i, dist);
				dists.push_back(temp);
			}
			
			sort(dists.begin(), dists.end(), compareSequenceDistance);  //sorts by distance lowest to highest
			
			//fill topmatches with numwanted closest sequences indexes
			for (int i = 0; i < numWanted; i++) {
				topMatches.push_back(dists[i].seq2);
			}
		
		}else{
			m->mothurOut("cannot find closest matches using distance method for " + query->getName() + " without aligned template sequences of the same length."); m->mothurOutEndLine();
			exit(1);
		}
		
		return topMatches;
	}
	catch(exception& e) {
		m->errorOut(e, "DistanceDB", "findClosestSequence");
		exit(1);
	}	
}
/**************************************************************************************************/
bool DistanceDB::isAligned(string seq){
	try {
		bool aligned;
		
		int pos = seq.find_first_of(".-");
		
		if (pos != seq.npos) {
			aligned = true;
		}else { aligned = false; }
		
		
		return aligned;
	}
	catch(exception& e) {
		m->errorOut(e, "DistanceDB", "isAligned");
		exit(1);
	}	
}

/**************************************************************************************************/
