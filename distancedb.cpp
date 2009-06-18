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


/**************************************************************************************************/

DistanceDB::DistanceDB(string fastaFileName, string distanceFileName) : Database(fastaFileName) {
	
	ifstream inputData;
	openInputFile(distanceFileName, inputData);
	int numCandSeqs=count(istreambuf_iterator<char>(inputData),istreambuf_iterator<char>(), '\n');	//	count the number of
	inputData.seekg(0);																		//	sequences

	hit closestMatch;
	string candidateSeqName;
	string junk;
	
	mostSimSequenceVector.resize(numCandSeqs);
	
	for(int i=0;i<numCandSeqs;i++){
		inputData >> candidateSeqName >> closestMatch.seqName >> closestMatch.indexNumber >> closestMatch.simScore;
//		getline(inputData, junk);	
		mostSimSequenceVector[i] = closestMatch;
	}
	cout << numCandSeqs << endl;
	searchIndex = 0;
	inputData.close();
}

/**************************************************************************************************/

Sequence DistanceDB::findClosestSequence(Sequence* candidateSeq){
	
	hit simAccession = mostSimSequenceVector[searchIndex];
//	string candidateSeqName, closestMatchSeqName, junk;
//	int closestMatchIndexNumber;
//	float closestMatchSimScore;
//	
//	inputData >> candidateSeqName >> closestMatchSeqName >> closestMatchIndexNumber >> closestMatchSimScore;
//	getline(inputData, junk);	

//	searchScore = 100. * closestMatchSimScore;

	searchScore = 100. * simAccession.simScore;
	searchIndex++;
//	return templateSequences[closestMatchIndexNumber];
	return templateSequences[simAccession.indexNumber];
	
}

/**************************************************************************************************/
