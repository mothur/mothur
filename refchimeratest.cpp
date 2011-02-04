/*
 *  refchimeratest.cpp
 *  Mothur
 *
 *  Created by Pat Schloss on 1/31/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "refchimeratest.h"
#include "mothur.h"

int MAXINT = numeric_limits<int>::max();

//***************************************************************************************************************

RefChimeraTest::RefChimeraTest(vector<Sequence>& refs, string chimeraReportFileName){

	m = MothurOut::getInstance();

	m->openOutputFile(chimeraReportFileName, chimeraReportFile);
	numRefSeqs = refs.size();

	referenceSeqs.resize(numRefSeqs);
	referenceNames.resize(numRefSeqs);
	for(int i=0;i<numRefSeqs;i++){
		referenceSeqs[i] = refs[i].getAligned();
		referenceNames[i] = refs[i].getName();
	}
	
	alignLength = referenceSeqs[0].length();

	chimeraReportFile << "queryName\tbestRef\tbestSequenceMismatch\tleftParentChi,rightParentChi\tbreakPointChi\tminMismatchToChimera\tdistToBestMera\tnumParents";
//	chimeraReportFile << "leftParentTri,middleParentTri,rightParentTri\tbreakPointTriA,breakPointTriB\tminMismatchToTrimera\tdistToBestMera\tnMera" << endl;

}

//***************************************************************************************************************

int RefChimeraTest::analyzeQuery(string queryName, string querySeq){
	
	vector<vector<int> > left(numRefSeqs);
	vector<int> singleLeft, bestLeft;
	vector<int> singleRight, bestRight;
	
	vector<vector<int> > right(numRefSeqs);
	for(int i=0;i<numRefSeqs;i++){
		left[i].assign(alignLength, 0);
	}
	right = left;
	
	int bestSequenceMismatch = getMismatches(querySeq, left, right, bestMatch);
	
	int leftParentBi, rightParentBi, breakPointBi;
	int minMismatchToChimera = getChimera(left, right, leftParentBi, rightParentBi, breakPointBi, singleLeft, bestLeft, singleRight, bestRight);
	
//	int minMismatchToTrimera = MAXINT;
//	int leftParentTri, middleParentTri, rightParentTri, breakPointTriA, breakPointTriB;
	
	int nMera = 0;
	string chimeraRefSeq = "";
	
	if(bestSequenceMismatch - minMismatchToChimera <= 3){
		nMera = 1;
		chimeraRefSeq = referenceSeqs[bestMatch];
	}
	else {

		nMera = 2;
		chimeraRefSeq = stitchBimera(leftParentBi, rightParentBi, breakPointBi);
		
//		minMismatchToTrimera = getTrimera(left, right, leftParentTri, middleParentTri, rightParentTri, breakPointTriA, breakPointTriB, singleLeft, bestLeft, singleRight, bestRight);
//
//		if(minMismatchToChimera - minMismatchToTrimera <= 3){
//			nMera = 2;
//			chimeraRefSeq = stitchBimera(leftParentBi, rightParentBi, breakPointBi);
//		}
//		else{			
//			nMera = 3;
//			chimeraRefSeq = stitchTrimera(leftParentTri, middleParentTri, rightParentTri, breakPointTriA, breakPointTriB);
//		}
		
	}
	double distToChimera = calcDistToChimera(querySeq, chimeraRefSeq);
	
//	double loonIndex = calcLoonIndex(querySeq, referenceSeqs[leftParentBi], referenceSeqs[rightParentBi], breakPointBi, binMatrix);		
	
	chimeraReportFile << queryName << '\t' << referenceNames[bestMatch] << '\t' << bestSequenceMismatch << '\t';
	chimeraReportFile << referenceNames[leftParentBi] << ',' << referenceNames[rightParentBi] << '\t' << breakPointBi << '\t';
	chimeraReportFile << minMismatchToChimera << '\t';
	
//	if(nMera == 1){
//		chimeraReportFile << "NA" << '\t' << "NA" << '\t' << "NA";
//	}
//	else{
//		chimeraReportFile << referenceNames[leftParentTri] << ',' << referenceNames[middleParentTri] << ',' << referenceNames[rightParentTri] << '\t' << breakPointTriA << ',' << breakPointTriB << '\t' << minMismatchToTrimera;	
//	}
	
	chimeraReportFile << '\t' << distToChimera << '\t' << nMera << endl;
		
	return nMera;
}

/**************************************************************************************************/

int RefChimeraTest::getMismatches(string& querySeq, vector<vector<int> >& left, vector<vector<int> >& right, int& bestRefSeq){
	
	int bestSequenceMismatch = MAXINT;
	
	for(int i=0;i<numRefSeqs;i++){
		
		int lDiffs = 0;
		for(int l=0;l<alignLength;l++){
			if(querySeq[l] != '.' && querySeq[l] != referenceSeqs[i][l]){
				lDiffs++;
			}
			left[i][l] = lDiffs;
		}
		
		int rDiffs = 0;
		int index = 0;
		for(int l=alignLength-1;l>=0;l--){
			if(querySeq[l] != '.' && querySeq[l] != referenceSeqs[i][l]){
				rDiffs++;
			}			
			right[i][index++] = rDiffs;
		}
		if(lDiffs < bestSequenceMismatch){
			bestSequenceMismatch = lDiffs;
			bestRefSeq = i;
		}
	}
	return bestSequenceMismatch;
}

/**************************************************************************************************/

int RefChimeraTest::getChimera(vector<vector<int> >& left, vector<vector<int> >& right, int& leftParent, int& rightParent, int& breakPoint, vector<int>& singleLeft, vector<int>& bestLeft, vector<int>& singleRight, vector<int>& bestRight){
	
	singleLeft.resize(alignLength, MAXINT);
	bestLeft.resize(alignLength, -1);
	
	for(int l=0;l<alignLength;l++){
		for(int i=0;i<numRefSeqs;i++){
			if(left[i][l] <= singleLeft[l]){
				singleLeft[l] = left[i][l];
				bestLeft[l] = i;
			}
		}
	}
	
	singleRight.resize(alignLength, MAXINT);
	bestRight.resize(alignLength, -1);
	
	for(int l=0;l<alignLength;l++){
		for(int i=0;i<numRefSeqs;i++){
			if(right[i][l] <= singleRight[l]){
				singleRight[l] = right[i][l];
				bestRight[l] = i;
			}
		}
	}
	
	int bestChimeraMismatches = MAXINT;
	leftParent = -1;
	rightParent = -1;
	breakPoint = -1;
	
	for(int l=0;l<alignLength;l++){
		int chimera = singleLeft[l] + singleRight[alignLength - l - 1];
		if(chimera < bestChimeraMismatches){
			bestChimeraMismatches = chimera;
			breakPoint = l;
			leftParent = bestLeft[l];
			rightParent = bestRight[alignLength - l - 1];
		}
	}
	
	return bestChimeraMismatches;
}

/**************************************************************************************************/

int RefChimeraTest::getTrimera(vector<vector<int> >& left, vector<vector<int> >& right, int& leftParent, int& middleParent, int& rightParent, int& breakPointA, int& breakPointB, vector<int>& singleLeft, vector<int>& bestLeft, vector<int>& singleRight, vector<int>& bestRight){
	
	int bestTrimeraMismatches = MAXINT;
	
	leftParent = -1;
	middleParent = -1;
	rightParent = -1;
	
	breakPointA = -1;
	breakPointB = -1;
	
	vector<vector<int> > minDelta(alignLength);
	vector<vector<int> > minDeltaSeq(alignLength);
	
	for(int i=0;i<alignLength;i++){
		minDelta[i].assign(alignLength, MAXINT);
		minDeltaSeq[i].assign(alignLength, -1);
	}
	
	for(int x=0;x<alignLength;x++){
		for(int y=x;y<alignLength-1;y++){
			
			for(int i=0;i<numRefSeqs;i++){
				int delta = left[i][y] - left[i][x];
				
				if(delta <= minDelta[x][y]){
					minDelta[x][y] = delta;
					minDeltaSeq[x][y] = i;					
				}				
			}
			minDelta[x][y] += singleLeft[x] + singleRight[alignLength - y - 2];
			
			if(minDelta[x][y] < bestTrimeraMismatches){
				bestTrimeraMismatches = minDelta[x][y];
				
				breakPointA = x;
				breakPointB = y;
				
				leftParent = bestLeft[x];
				middleParent = minDeltaSeq[x][y];
				rightParent = bestRight[alignLength - y - 2];				
			}
		}		
	}
	return bestTrimeraMismatches;
}

/**************************************************************************************************/

string RefChimeraTest::stitchBimera(int leftParent, int rightParent, int breakPoint){
	
	string chimeraRefSeq = referenceSeqs[leftParent].substr(0, breakPoint) + referenceSeqs[rightParent].substr(breakPoint);
	return chimeraRefSeq;
	
}

/**************************************************************************************************/

string RefChimeraTest::stitchTrimera(int leftParent, int middleParent, int rightParent, int breakPointA, int breakPointB){
	
	string chimeraRefSeq = referenceSeqs[leftParent].substr(0, breakPointA) + referenceSeqs[middleParent].substr(breakPointA, breakPointB-breakPointA) + referenceSeqs[rightParent].substr(breakPointB);
	
	return chimeraRefSeq;
}

/**************************************************************************************************/

double RefChimeraTest::calcDistToChimera(string& querySeq, string& chimeraRefSeq){
	
	int match = 0;
	int mismatch = 0;
	
	for(int i=0;i<alignLength;i++){
		if(querySeq[i] != '.' && chimeraRefSeq[i] != '.'){
			if(querySeq[i] == '-' && chimeraRefSeq[i] == '-'){	/*	do nothing	*/	}
			else if(querySeq[i] == chimeraRefSeq[i]){
				match++;
			}
			else{
				mismatch++;
			}			
		}
	}
	
	return (double)mismatch / (double)(mismatch + match);	
}

//***************************************************************************************************************

int RefChimeraTest::getClosestRefIndex(){

	return bestMatch;
	
}

//***************************************************************************************************************
