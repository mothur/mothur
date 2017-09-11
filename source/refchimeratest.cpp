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

RefChimeraTest::RefChimeraTest(vector<Sequence>& refs, bool aligned) : aligned(aligned){

	m = MothurOut::getInstance();

	numRefSeqs = refs.size();

	referenceSeqs.resize(numRefSeqs);
	referenceNames.resize(numRefSeqs);
	for(int i=0;i<numRefSeqs;i++){
		if (aligned) { referenceSeqs[i] = refs[i].getAligned(); }
        else { referenceSeqs[i] = refs[i].getUnaligned(); }
		referenceNames[i] = refs[i].getName();
	}
	
	alignLength = referenceSeqs[0].length();
    bestMatch = 0;

}
//***************************************************************************************************************

int RefChimeraTest::printHeader(ofstream& chimeraReportFile){
	try {
		chimeraReportFile << "queryName\tbestRef\tbestSequenceMismatch\tleftParentChi,rightParentChi\tbreakPointChi\tminMismatchToChimera\tdistToBestMera\tnumParents" << endl;
		return 0; 
	}catch(exception& e) {
		m->errorOut(e, "RefChimeraTest", "printHeader");
		exit(1);
	}
}

//***************************************************************************************************************

int RefChimeraTest::analyzeQuery(string queryName, string querySeq, ofstream& chimeraReportFile){

    int numParents = -1;
    
    if(aligned){
        numParents = analyzeAlignedQuery(queryName, querySeq, chimeraReportFile);
    }
    else{
        numParents = analyzeUnalignedQuery(queryName, querySeq, chimeraReportFile);
    }
    
    return numParents;
    
}
 
//***************************************************************************************************************

int RefChimeraTest::analyzeAlignedQuery(string queryName, string querySeq, ofstream& chimeraReportFile){
	
    vector<vector<int> > left; left.resize(numRefSeqs);
    vector<vector<int> > right; right.resize(numRefSeqs);
	vector<int> singleLeft, bestLeft;
	vector<int> singleRight, bestRight;
	
	for(int i=0;i<numRefSeqs;i++){
		left[i].assign(alignLength, 0);
        right[i].assign(alignLength, 0);
	}
	
    int bestMatchIndex;
	int bestSequenceMismatch = getAlignedMismatches(querySeq, left, right, bestMatchIndex);
	
	int leftParentBi, rightParentBi, breakPointBi;
	int minMismatchToChimera = getChimera(left, right, leftParentBi, rightParentBi, breakPointBi, singleLeft, bestLeft, singleRight, bestRight);
	
	int nMera = 0;
	string chimeraRefSeq = "";
	
	if(bestSequenceMismatch - minMismatchToChimera >= 3){// || (minMismatchToChimera == 0 && bestSequenceMismatch != 0)){
		nMera = 2;
		chimeraRefSeq = stitchBimera(leftParentBi, rightParentBi, breakPointBi);
	}
	else{
		nMera = 1;
		chimeraRefSeq = referenceSeqs[bestMatchIndex];
	}
	
	bestRefAlignment = chimeraRefSeq;
    bestQueryAlignment = querySeq;
    
	double distToChimera = calcDistToChimera(bestQueryAlignment, bestRefAlignment);
	
	chimeraReportFile << queryName << '\t' << referenceNames[bestMatchIndex] << '\t' << bestSequenceMismatch << '\t';
	chimeraReportFile << referenceNames[leftParentBi] << ',' << referenceNames[rightParentBi] << '\t' << breakPointBi << '\t';
	chimeraReportFile << minMismatchToChimera << '\t';
    chimeraReportFile << '\t' << distToChimera << '\t' << nMera << endl;
    
    bestMatch = bestMatchIndex;
		
	return nMera;
}

//***************************************************************************************************************

int RefChimeraTest::analyzeUnalignedQuery(string queryName, string querySeq, ofstream& chimeraReportFile){
	
    int nMera = 0;
    
    int seqLength = querySeq.length();
    
    vector<string> queryAlign; queryAlign.resize(numRefSeqs);
    vector<string> refAlign; refAlign.resize(numRefSeqs);

    vector<vector<int> > leftDiffs; leftDiffs.resize(numRefSeqs);
    vector<vector<int> > rightDiffs; rightDiffs.resize(numRefSeqs);
    vector<vector<int> > leftMaps; leftMaps.resize(numRefSeqs);
    vector<vector<int> > rightMaps; rightMaps.resize(numRefSeqs);
    
    int bestRefIndex = -1;
    int bestRefDiffs = numeric_limits<int>::max();
    double bestRefLength = 0;
    
    for(int i=0;i<numRefSeqs;i++){
        double length = 0;
        double diffs = alignQueryToReferences(querySeq, referenceSeqs[i], queryAlign[i], refAlign[i], length);
        if(diffs < bestRefDiffs){
            bestRefDiffs = diffs;
            bestRefLength = length;
            bestRefIndex = i;
        }
    }

    if(bestRefDiffs >= 3){
        for(int i=0;i<numRefSeqs;i++){
            leftDiffs[i].assign(seqLength, 0);
            rightDiffs[i].assign(seqLength, 0);
            leftMaps[i].assign(seqLength, 0);
            rightMaps[i].assign(seqLength, 0);
            
            getUnalignedDiffs(queryAlign[i], refAlign[i], leftDiffs[i], leftMaps[i], rightDiffs[i], rightMaps[i]);
        }
    
        vector<int> singleLeft(seqLength, numeric_limits<int>::max());
        vector<int> bestLeft(seqLength, -1);
        
        for(int l=0;l<seqLength;l++){
            
            for(int i=0;i<numRefSeqs;i++){            
                if(leftDiffs[i][l] < singleLeft[l]){
                    singleLeft[l] = leftDiffs[i][l];
                    bestLeft[l] = i;
                }
            }
        }
        
        vector<int> singleRight(seqLength, numeric_limits<int>::max());
        vector<int> bestRight(seqLength, -1);
        
        for(int l=0;l<seqLength;l++){
                
            for(int i=0;i<numRefSeqs;i++){
                if(rightDiffs[i][l] < singleRight[l]){
                    singleRight[l] = rightDiffs[i][l];
                    bestRight[l] = i;
                }
            }
        }
        
        int bestChimeraMismatches = numeric_limits<int>::max();
        int leftParent = 0;
        int rightParent = 0;
        int breakPoint = 0;
        
        for(int l=0;l<seqLength-1;l++){
            
            int chimera = singleLeft[l] + singleRight[seqLength - l - 2];
            if(chimera < bestChimeraMismatches){
                bestChimeraMismatches = chimera;
                breakPoint = l;
                leftParent = bestLeft[l];
                rightParent = bestRight[seqLength - l - 2];
            }
        }
        
        string reference;
        
        if(bestRefDiffs - bestChimeraMismatches >= 3){// || (minMismatchToChimera == 0 && bestSequenceMismatch != 0)){
            nMera = 2;
            
            int breakLeft = leftMaps[leftParent][breakPoint];
            int breakRight = rightMaps[rightParent][rightMaps[rightParent].size() - breakPoint - 2];
            
            string left = refAlign[leftParent];
            string right = refAlign[rightParent];
            
            for(int i=0;i<=breakLeft;i++){
                
                if (m->getControl_pressed()) { return 0; }
                
                if(left[i] != '-' && left[i] != '.'){
                    reference += left[i];
                }
            }
            
            
            for(int i=breakRight;i<right.length();i++){
                
                if (m->getControl_pressed()) { return 0; }
                
                if(right[i] != '-' && right[i] != '.'){
                    reference += right[i];
                }
            }

        }
        else{
            nMera = 1;
            reference = referenceSeqs[bestRefIndex];
        }

        double alignLength;
        double finalDiffs = alignQueryToReferences(querySeq, reference, bestQueryAlignment, bestRefAlignment, alignLength);
        double finalDistance = finalDiffs / alignLength;

        chimeraReportFile << queryName << '\t' << referenceNames[bestRefIndex] << '\t' << bestRefDiffs << '\t';
        chimeraReportFile << referenceNames[leftParent] << ',' << referenceNames[rightParent] << '\t' << breakPoint << '\t';
        chimeraReportFile << bestChimeraMismatches << '\t';
        chimeraReportFile << '\t' << finalDistance << '\t' << nMera << endl;
    }
    else{
        bestQueryAlignment = queryAlign[bestRefIndex];
        bestRefAlignment = refAlign[bestRefIndex];
        nMera = 1;
        
        chimeraReportFile << queryName << '\t' << referenceNames[bestRefIndex] << '\t' << bestRefDiffs << '\t';
        chimeraReportFile << "NA\tNA\tNA\tNA\t1" << endl;
    } 
    
    bestMatch = bestRefIndex;
    return nMera;
}

/**************************************************************************************************/

double RefChimeraTest::alignQueryToReferences(string query, string reference, string& qAlign, string& rAlign, double& length){
    
    
    try {
		double GAP = -5;
		double MATCH = 1;
		double MISMATCH = -1;
		
		int queryLength = query.length();
		int refLength = reference.length();
		
        vector<vector<double> > alignMatrix; alignMatrix.resize(queryLength + 1);
        vector<vector<char> > alignMoves; alignMoves.resize(queryLength + 1);
		
		for(int i=0;i<=queryLength;i++){
			if (m->getControl_pressed()) { return 0; }
			alignMatrix[i].resize(refLength + 1, 0);
			alignMoves[i].resize(refLength + 1, 'x');
		}
		
		for(int i=0;i<=queryLength;i++){
			if (m->getControl_pressed()) { return 0; }
			alignMatrix[i][0] = 0;//GAP * i;
			alignMoves[i][0] = 'u';
		}
		
		for(int i=0;i<=refLength;i++){
			if (m->getControl_pressed()) { return 0; }
			alignMatrix[0][i] = 0;//GAP * i;
			alignMoves[0][i] = 'l';
		}
		        
		for(int i=1;i<=queryLength;i++){
			
			if (m->getControl_pressed()) { return 0; }
			
			for(int j=1;j<=refLength;j++){
				
				double nogapScore;		
				if(query[i-1] == reference[j-1]){	nogapScore = alignMatrix[i-1][j-1] + MATCH;		}
				else							{	nogapScore = alignMatrix[i-1][j-1] + MISMATCH;	}
				
				double leftScore;
				if(i == queryLength)			{	leftScore = alignMatrix[i][j-1];				}
				else							{	leftScore = alignMatrix[i][j-1] + GAP;			}
				
				
				double upScore;
				if(j == refLength)				{	upScore = alignMatrix[i-1][j];					}
				else							{	upScore = alignMatrix[i-1][j] + GAP;			}
				
				if(nogapScore > leftScore){
					if(nogapScore > upScore){
						alignMoves[i][j] = 'd';
						alignMatrix[i][j] = nogapScore;
					}
					else{
						alignMoves[i][j] = 'u';
						alignMatrix[i][j] = upScore;
					}
				}
				else{
					if(leftScore > upScore){
						alignMoves[i][j] = 'l';
						alignMatrix[i][j] = leftScore;
					}
					else{
						alignMoves[i][j] = 'u';
						alignMatrix[i][j] = upScore;
					}
				}
			}
		}
        		
		int end = refLength - 1;
        int maxRow = 0;
        double maxRowValue = -2147483647;
        for(int i=0;i<queryLength;i++){
            if(alignMatrix[i][end] > maxRowValue){
                maxRow = i;
                maxRowValue = alignMatrix[i][end];
            }
        }
        
        end = queryLength - 1;
        int maxColumn = 0;
        double maxColumnValue = -2147483647;

        for(int j=0;j<refLength;j++){
            if(alignMatrix[end][j] > maxColumnValue){
                maxColumn = j;
                maxColumnValue = alignMatrix[end][j];
            }
        }

        int row = queryLength-1;
        int column = refLength-1;
        
        if(maxColumn == column && maxRow == row){}	//	if the max values are the lower right corner, then we're good
        else if(alignMatrix[row][maxColumn] < alignMatrix[maxRow][column]){
            for(int i=maxRow+1;i<queryLength;i++){			//	decide whether sequence A or B needs the gaps at the end either set 
                alignMoves[i][column] = 'u';//	the pointer upwards or...
            }
            
        }
        else {
            for(int i=maxColumn+1;i<refLength;i++){
                alignMoves[row][i] = 'l';	//	...to the left
            }
        }

        int i = queryLength;
		int j = refLength;

        
		qAlign = "";
		rAlign = "";
        
		int diffs = 0;
		length = 0;

		while(i > 0 && j > 0){
			
			if (m->getControl_pressed()) { return 0; }
			
			if(alignMoves[i][j] == 'd'){
				qAlign = query[i-1] + qAlign;
				rAlign = reference[j-1] + rAlign;
                
				if(query[i-1] != reference[j-1]){	diffs++;	}
				length++;
				
				i--;
				j--;
			}
			else if(alignMoves[i][j] == 'u'){
				qAlign = query[i-1] + qAlign;
				
				if(j != refLength)	{	rAlign = '-' + rAlign;	diffs++;	length++;	}
				else				{	rAlign = '.' + rAlign;	}
				i--;
			}
			else if(alignMoves[i][j] == 'l'){
				rAlign = reference[j-1] + rAlign;
				
				if(i != queryLength){	qAlign = '-' + qAlign;	diffs++;	length++;	}
				else				{	qAlign = '.' + qAlign;	}
				j--;
			}
		}

        if(i>0){
            qAlign = query.substr(0, i) + qAlign;
            rAlign = string(i, '.') + rAlign;
        }
		else if(j>0){
            qAlign = string(j, '.') + qAlign;
            rAlign = reference.substr(0, j) + rAlign;
        }

        
		return diffs;
	}
	catch(exception& e) {
		m->errorOut(e, "RefChimeraTest", "alignQueryToReferences");
		exit(1);
	}
}

/**************************************************************************************************/

int RefChimeraTest::getUnalignedDiffs(string qAlign, string rAlign, vector<int>& leftDiffs, vector<int>& leftMap, vector<int>& rightDiffs, vector<int>& rightMap){
	try {
		int alignLength = qAlign.length();
		
		int lDiffs = 0;
		int lCount = 0;
		for(int l=0;l<alignLength;l++){
			
			if (m->getControl_pressed()) { return 0; }
			
			if(qAlign[l] == '-'){
				lDiffs++;		
			}
			else if(qAlign[l] != '.'){
				
				if(rAlign[l] == '-'){
					lDiffs++;
				}
				else if(qAlign[l] != rAlign[l]){;// && rAlign[l] != '.'){
					lDiffs++;
				}
				leftDiffs[lCount] = lDiffs;
				leftMap[lCount] = l;
				
				lCount++;
			}			
		}

		int rDiffs = 0;
		int rCount = 0;
		for(int l=alignLength-1;l>=0;l--){
			
			if (m->getControl_pressed()) { return 0; }
			
			if(qAlign[l] == '-'){
				rDiffs++;		
			}
			else if(qAlign[l] != '.'){
				
				if(rAlign[l] == '-'){
					rDiffs++;
				}
				else if(qAlign[l] != rAlign[l]){;// && rAlign[l] != '.'){
					rDiffs++;
				}
				
				rightDiffs[rCount] = rDiffs;
				rightMap[rCount] = l;
				rCount++;
			}
			
		}

		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "RefChimeraTest", "getUnalignedDiffs");
		exit(1);
	}
}

/**************************************************************************************************/

int RefChimeraTest::getAlignedMismatches(string& querySeq, vector<vector<int> >& left, vector<vector<int> >& right, int& bestRefSeq){
	
	int bestSequenceMismatch = MAXINT;
	
	for(int i=0;i<numRefSeqs;i++){
		
		int lDiffs = 0;
		
		for(int l=0;l<alignLength;l++){
			if(querySeq[l] != '.' && referenceSeqs[i][l] != '.' && querySeq[l] != referenceSeqs[i][l] && referenceSeqs[i][l] != 'N'){
				lDiffs++;
			}
			left[i][l] = lDiffs;
		}
		
		int rDiffs = 0;
		int index = 0;
		for(int l=alignLength-1;l>=0;l--){
			if(querySeq[l] != '.' && referenceSeqs[i][l] != '.' && querySeq[l] != referenceSeqs[i][l] && referenceSeqs[i][l] != 'N'){
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
	
    vector<vector<int> > minDelta; minDelta.resize(alignLength);
    vector<vector<int> > minDeltaSeq; minDeltaSeq.resize(alignLength);
	
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
//		if(querySeq[i] != '.' && chimeraRefSeq[i] != '.'){
		if(chimeraRefSeq[i] != '.' && querySeq[i] != '.'){
			if(querySeq[i] == '-' && chimeraRefSeq[i] == '-' && chimeraRefSeq[i] != 'N'){	/*	do nothing	*/	}
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

string RefChimeraTest::getQueryAlignment(){
    
	return bestQueryAlignment;
	
}

//***************************************************************************************************************

string RefChimeraTest::getClosestRefAlignment(){
    
	return bestRefAlignment;
	
}

//***************************************************************************************************************
