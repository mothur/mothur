/*
 *  myPerseus.cpp
 *  
 *
 *  Created by Pat Schloss on 9/5/11.
 *  Copyright 2011 Patrick D. Schloss. All rights reserved.
 *
 */

#include "myPerseus.h"

/**************************************************************************************************/
int PERSEUSMAXINT = numeric_limits<int>::max();
/**************************************************************************************************/

vector<vector<double> > Perseus::binomial(int maxOrder){
	try {
		vector<vector<double> > binomial(maxOrder+1);
		
		for(int i=0;i<=maxOrder;i++){
			binomial[i].resize(maxOrder+1);
			binomial[i][0]=1;
			binomial[0][i]=0;
		}
		binomial[0][0]=1;
		
		binomial[1][0]=1;
		binomial[1][1]=1;
		
		for(int i=2;i<=maxOrder;i++){
			binomial[1][i]=0;
		}
		
		for(int i=2;i<=maxOrder;i++){
			for(int j=1;j<=maxOrder;j++){
				if(i==j){	binomial[i][j]=1;									}
				if(j>i)	{	binomial[i][j]=0;									}
				else	{	binomial[i][j]=binomial[i-1][j-1]+binomial[i-1][j];	}
			}
		}
		
		return binomial;
	}
	catch(exception& e) {
		m->errorOut(e, "Perseus", "binomial");
		exit(1);
	}
}

/**************************************************************************************************/
double Perseus::basicPairwiseAlignSeqs(string query, string reference, string& qAlign, string& rAlign, pwModel model){
	try {
		double GAP = model.GAP_OPEN;
		double MATCH = model.MATCH;
		double MISMATCH = model.MISMATCH;
		
		int queryLength = query.size();
		int refLength = reference.size();
		
		vector<vector<double> > alignMatrix(queryLength + 1);
		vector<vector<char> > alignMoves(queryLength + 1);
		
		for(int i=0;i<=queryLength;i++){
			if (m->getControl_pressed()) { return 0; }
			alignMatrix[i].resize(refLength + 1, 0);
			alignMoves[i].resize(refLength + 1, 'x');
		}
		
		for(int i=0;i<=queryLength;i++){
			if (m->getControl_pressed()) { return 0; }
			alignMatrix[i][0] = GAP * i;
			alignMoves[i][0] = 'u';
		}
		
		for(int i=0;i<=refLength;i++){
			if (m->getControl_pressed()) { return 0; }
			alignMatrix[0][i] = GAP * i;
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
		
		int i = queryLength;
		int j = refLength;
		
		qAlign = "";
		rAlign = "";
			
		int diffs = 0;
		int length = 0;
		
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
		
		while(i>0){
			
			if (m->getControl_pressed()) { return 0; }
			
			rAlign = '.' + rAlign;
			qAlign = query[i-1] + qAlign;
			i--;
		}
		
		while(j>0){
			
			if (m->getControl_pressed()) { return 0; }
			
			rAlign = reference[j-1] + rAlign;
			qAlign = '.' + qAlign;
			j--;
		}
		
		

		return double(diffs)/double(length);
	}
	catch(exception& e) {
		m->errorOut(e, "Perseus", "basicPairwiseAlignSeqs");
		exit(1);
	}
	
}
/**************************************************************************************************/
int Perseus::getDiffs(string qAlign, string rAlign, vector<int>& leftDiffs, vector<int>& leftMap, vector<int>& rightDiffs, vector<int>& rightMap){
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
				else if(qAlign[l] != rAlign[l] && rAlign[l] != '.'){
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
				else if(qAlign[l] != rAlign[l] && rAlign[l] != '.'){
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
		m->errorOut(e, "Perseus", "getDiffs");
		exit(1);
	}
}
/**************************************************************************************************/
int Perseus::getLastMatch(char direction, vector<vector<char> >& alignMoves, int i, int j, string& seqA, string& seqB){
	try {
		char nullReturn = -1;
		
		while(i>=1 && j>=1){
			
			if (m->getControl_pressed()) { return 0; }
			
			if(direction == 'd'){
				if(seqA[i-1] == seqB[j-1])	{	return seqA[i-1];	}
				else						{	return nullReturn;	}
			}
			
			else if(direction == 'l')		{	j--;				}
			else							{	i--;				}
			
			direction = alignMoves[i][j];
		}
		
		return nullReturn;
	}
	catch(exception& e) {
		m->errorOut(e, "Perseus", "getLastMatch");
		exit(1);
	}
}

/**************************************************************************************************/

int Perseus::toInt(char b){
	try {
		if(b == 'A')		{	return 0;	}
		else if(b == 'C')	{	return 1;	}
		else if(b == 'T')	{	return 2;	}
		else if(b == 'G')	{	return 3;	}
		else { m->mothurOut("[ERROR]: " + toString(b) + " is not ATGC.\n");  return -1; }
	}
	catch(exception& e) {
		m->errorOut(e, "Perseus", "toInt");
		exit(1);
	}
}

/**************************************************************************************************/

double Perseus::modeledPairwiseAlignSeqs(string query, string reference, string& qAlign, string& rAlign, vector<vector<double> >& correctMatrix){
	try {
		int queryLength = query.size();
		int refLength = reference.size();
		
		vector<vector<double> > alignMatrix(queryLength + 1);
		vector<vector<char> > alignMoves(queryLength + 1);
		
		for(int i=0;i<=queryLength;i++){
			if (m->getControl_pressed()) { return 0; }
			alignMatrix[i].resize(refLength + 1, 0);
			alignMoves[i].resize(refLength + 1, 'x');
		}
		
		for(int i=0;i<=queryLength;i++){
			if (m->getControl_pressed()) { return 0; }
			alignMatrix[i][0] = 15.0 * i;
			alignMoves[i][0] = 'u';
		}
		
		for(int i=0;i<=refLength;i++){
			if (m->getControl_pressed()) { return 0; }
			alignMatrix[0][i] = 15.0 * i;
			alignMoves[0][i] = 'l';
		}
		
		for(int i=1;i<=queryLength;i++){
			
			if (m->getControl_pressed()) { return 0; }
			
			for(int j=1;j<=refLength;j++){
				
				double nogap;		
				nogap = alignMatrix[i-1][j-1] + correctMatrix[toInt(query[i-1])][toInt(reference[j-1])];			
				
				double gap;
				
				double left;
				if(i == queryLength){ //terminal gap
					left = alignMatrix[i][j-1];
				}
				else{
					if(reference[j-1] == getLastMatch('l', alignMoves, i, j, query, reference)){
						gap = 4.0;
					}
					else{
						gap = 15.0;
					}
					
					left = alignMatrix[i][j-1] + gap;
				}
				
				double up;
				if(j == refLength){ //terminal gap
					up = alignMatrix[i-1][j];
				}
				else{
					
					if(query[i-1] == getLastMatch('u', alignMoves, i, j, query, reference)){
						gap = 4.0;
					}
					else{
						gap = 15.0;
					}
					
					up = alignMatrix[i-1][j] + gap;
				}
				
				
				if(nogap < left){
					if(nogap < up){
						alignMoves[i][j] = 'd';
						alignMatrix[i][j] = nogap;
					}
					else{
						alignMoves[i][j] = 'u';
						alignMatrix[i][j] = up;
					}
				}
				else{
					if(left < up){
						alignMoves[i][j] = 'l';
						alignMatrix[i][j] = left;
					}
					else{
						alignMoves[i][j] = 'u';
						alignMatrix[i][j] = up;
					}
				}
			}
		}

		int i = queryLength;
		int j = refLength;
		
		int alignLength = 0;
		
		while(i > 0 && j > 0){
			
			if (m->getControl_pressed()) { return 0; }
			
			if(alignMoves[i][j] == 'd'){
				qAlign = query[i-1] + qAlign;
				rAlign = reference[j-1] + rAlign;
				alignLength++;
				i--;
				j--;
			}
			else if(alignMoves[i][j] == 'u'){
				if(j != refLength){
					qAlign = query[i-1] + qAlign;
					rAlign = '-' + rAlign;
					alignLength++;
				}
				
				i--;
			}
			else if(alignMoves[i][j] == 'l'){
				if(i != queryLength){
					qAlign = '-' + qAlign;
					rAlign = reference[j-1] + rAlign;
					alignLength++;				
				}
				
				j--;
			}
		}

		return alignMatrix[queryLength][refLength] / (double)alignLength;
	}
	catch(exception& e) {
		m->errorOut(e, "Perseus", "modeledPairwiseAlignSeqs");
		exit(1);
	}
}

/**************************************************************************************************/
int Perseus::getAlignments(int curSequenceIndex, vector<seqData> sequences, vector<pwAlign>& alignments, vector<vector<int> >& leftDiffs, vector<vector<int> >& leftMaps, vector<vector<int> >& rightDiffs, vector<vector<int> >& rightMaps, int& bestRefSeq, int& bestRefDiff, vector<bool>& restricted){
	try {
		int numSeqs = sequences.size();
		//int bestSequenceMismatch = PERSEUSMAXINT;

		string curSequence = sequences[curSequenceIndex].sequence;
		int curFrequency = sequences[curSequenceIndex].frequency; 

		bestRefSeq = -1;
		
		int bestIndex = -1;
		int bestDiffs = PERSEUSMAXINT;
		int comparisons = 0;
			
		pwModel model(0, -1, -1.5);
		
		for(int i=0;i<numSeqs;i++){
			
			if (m->getControl_pressed()) { return 0; }
			
			if(i != curSequenceIndex && restricted[i] != 1 && sequences[i].frequency >= 2 * curFrequency){
				string refSequence = sequences[i].sequence;
				
				leftDiffs[i].assign(curSequence.length(), 0);
				leftMaps[i].assign(curSequence.length(), 0);
				rightDiffs[i].assign(curSequence.length(), 0);
				rightMaps[i].assign(curSequence.length(), 0);
				
				basicPairwiseAlignSeqs(curSequence, refSequence, alignments[i].query, alignments[i].reference, model);
				

				getDiffs(alignments[i].query, alignments[i].reference, leftDiffs[i], leftMaps[i], rightDiffs[i], rightMaps[i]);
				
				int diffs = rightDiffs[i][curSequence.length()-1];							

				if(diffs < bestDiffs){
					bestDiffs = diffs;
					bestIndex = i;
				}
				comparisons++;
				restricted[i] = 0;
			}
			else{
				restricted[i] = 1;
			}
		}

		bestRefSeq = bestIndex;
		bestRefDiff = bestDiffs;
		
		return comparisons;
	}
	catch(exception& e) {
		m->errorOut(e, "Perseus", "getAlignments");
		exit(1);
	}
}
/**************************************************************************************************/
int Perseus::getChimera(vector<seqData> sequences,
			   vector<vector<int> >& leftDiffs, 
			   vector<vector<int> >& rightDiffs,
			   int& leftParent, 
			   int& rightParent, 
			   int& breakPoint,
			   vector<int>& singleLeft, 
			   vector<int>& bestLeft, 
			   vector<int>& singleRight, 
			   vector<int>& bestRight, 
			   vector<bool> restricted){
	try {
		int numRefSeqs = restricted.size();
		int seqLength = leftDiffs[0].size();
		
		singleLeft.resize(seqLength, PERSEUSMAXINT);
		bestLeft.resize(seqLength, -1);

		for(int l=0;l<seqLength;l++){
			
			if (m->getControl_pressed()) { return 0; }
			
			for(int i=0;i<numRefSeqs;i++){
				
				if(!restricted[i]){
					if(((leftDiffs[i][l] < singleLeft[l]) && sequences[i].frequency) || ((leftDiffs[i][l] == singleLeft[l]) && (sequences[i].frequency > sequences[bestLeft[l]].frequency))){
						singleLeft[l] = leftDiffs[i][l];
						bestLeft[l] = i;
					}
				}
			}
		}
		
		singleRight.resize(seqLength, PERSEUSMAXINT);
		bestRight.resize(seqLength, -1);
		
		for(int l=0;l<seqLength;l++){
			
			if (m->getControl_pressed()) { return 0; }
			
			for(int i=0;i<numRefSeqs;i++){
				
				if(!restricted[i]){
					if((rightDiffs[i][l] < singleRight[l] && sequences[i].frequency) || ((rightDiffs[i][l] == singleRight[l] && sequences[i].frequency > sequences[bestRight[l]].frequency))){
						singleRight[l] = rightDiffs[i][l];
						bestRight[l] = i;
					}
				}
			}
		}

		
		
		int bestChimeraMismatches = PERSEUSMAXINT;
		leftParent = -1;
		rightParent = -1;
		breakPoint = -1;
		
		for(int l=0;l<seqLength-1;l++){
			
			if (m->getControl_pressed()) { return 0; }
			
			int chimera = singleLeft[l] + singleRight[seqLength - l - 2];
			if(chimera < bestChimeraMismatches){
				bestChimeraMismatches = chimera;
				breakPoint = l;
				leftParent = bestLeft[l];
				rightParent = bestRight[seqLength - l - 2];
			}
		}
		
		return bestChimeraMismatches;
	}
	catch(exception& e) {
		m->errorOut(e, "Perseus", "getChimera");
		exit(1);
	}
}

/**************************************************************************************************/

string Perseus::stitchBimera(vector<pwAlign>& alignments, int leftParent, int rightParent, int breakPoint, vector<vector<int> >& leftMaps, vector<vector<int> >& rightMaps){
	try {
		int breakLeft = leftMaps[leftParent][breakPoint];
		int breakRight = rightMaps[rightParent][rightMaps[rightParent].size() - breakPoint - 2];
		
		string left = alignments[leftParent].reference;
		string right = alignments[rightParent].reference;
		string chimera = "";
		
		for(int i=0;i<=breakLeft;i++){
			
			if (m->getControl_pressed()) { return 0; }
			
			if(left[i] != '-' && left[i] != '.'){
				chimera += left[i];
			}
		}
		

		for(int i=breakRight;i<right.length();i++){
			
			if (m->getControl_pressed()) { return 0; }
			
			if(right[i] != '-' && right[i] != '.'){
				chimera += right[i];
			}
		}
		
		return chimera;
	}
	catch(exception& e) {
		m->errorOut(e, "Perseus", "stitchBimera");
		exit(1);
	}
}
/**************************************************************************************************/
int Perseus::getTrimera(vector<seqData>& sequences,
			   vector<vector<int> >& leftDiffs,
			   int& leftParent,
			   int& middleParent,
			   int& rightParent,
			   int& breakPointA,
			   int& breakPointB,
			   vector<int>& singleLeft,
			   vector<int>& bestLeft, 
			   vector<int>& singleRight,
			   vector<int>& bestRight,
			   vector<bool> restricted){
	try {
		int numRefSeqs = leftDiffs.size();
		int alignLength = leftDiffs[0].size();
		int bestTrimeraMismatches = PERSEUSMAXINT;
		
		leftParent = -1;
		middleParent = -1;
		rightParent = -1;
		
		breakPointA = -1;
		breakPointB = -1;
		
		vector<vector<int> > minDelta(alignLength);
		vector<vector<int> > minDeltaSeq(alignLength);
		
		for(int i=0;i<alignLength;i++){
			if (m->getControl_pressed()) { return 0; }
			minDelta[i].assign(alignLength, PERSEUSMAXINT);
			minDeltaSeq[i].assign(alignLength, -1);
		}
		
		for(int x=0;x<alignLength;x++){
			for(int y=x;y<alignLength-1;y++){
				for(int i=0;i<numRefSeqs;i++){
					
					if (m->getControl_pressed()) { return 0; }
					
					if(!restricted[i]){
						int delta = leftDiffs[i][y] - leftDiffs[i][x];

						if(delta < minDelta[x][y] || (delta == minDelta[x][y] && sequences[i].frequency > sequences[minDeltaSeq[x][y]].frequency)){
							minDelta[x][y] = delta;
							minDeltaSeq[x][y] = i;					
						}				
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
	catch(exception& e) {
		m->errorOut(e, "Perseus", "getTrimera");
		exit(1);
	}
}

/**************************************************************************************************/

string Perseus::stitchTrimera(vector<pwAlign> alignments, int leftParent, int middleParent, int rightParent, int breakPointA, int breakPointB, vector<vector<int> >& leftMaps, vector<vector<int> >& rightMaps){
	try {
		int p1SplitPoint = leftMaps[leftParent][breakPointA];
		int p2SplitPoint = leftMaps[middleParent][breakPointB];
		int p3SplitPoint = rightMaps[rightParent][rightMaps[rightParent].size() - breakPointB - 2];
		
		string chimeraRefSeq;
		for(int i=0;i<=p1SplitPoint;i++){
			if (m->getControl_pressed()) { return chimeraRefSeq; }
			if(alignments[leftParent].reference[i] != '-' && alignments[leftParent].reference[i] != '.'){
				chimeraRefSeq += alignments[leftParent].reference[i];
			}
		}
		
		for(int i=p1SplitPoint+1;i<=p2SplitPoint;i++){
			if (m->getControl_pressed()) { return chimeraRefSeq; }
			if(alignments[middleParent].reference[i] != '-' && alignments[middleParent].reference[i] != '.'){
				chimeraRefSeq += alignments[middleParent].reference[i];
			}
		}
		
		for(int i=p3SplitPoint;i<alignments[rightParent].reference.length();i++){
			if (m->getControl_pressed()) { return chimeraRefSeq; }
			if(alignments[rightParent].reference[i] != '-' && alignments[rightParent].reference[i] != '.'){
				chimeraRefSeq += alignments[rightParent].reference[i];
			}
		}

		return chimeraRefSeq;
	}
	catch(exception& e) {
		m->errorOut(e, "Perseus", "stitchTrimera");
		exit(1);
	}
}

/**************************************************************************************************/

int Perseus::threeWayAlign(string query, string parent1, string parent2, string& qAlign, string& aAlign, string& bAlign){
	try {
		pwModel model(1.0, -1.0, -5.0);
		
		string qL, rL;
		string qR, rR;

		basicPairwiseAlignSeqs(query, parent1, qL, rL, model);	
		basicPairwiseAlignSeqs(query, parent2, qR, rR, model);

		int lLength = qL.length();
		int rLength = qR.length();
		
		string qLNew, rLNew;
		string qRNew, rRNew;
		
		int lIndex = 0;
		int rIndex = 0;
		
		while(lIndex<lLength || rIndex<rLength){
			
			if (m->getControl_pressed()) { return 0; }
			
			if(qL[lIndex] == qR[rIndex]){
				qLNew += qL[lIndex];
				rLNew += rL[lIndex];
				lIndex++;

				qRNew += qR[rIndex];
				rRNew += rR[rIndex];
				rIndex++;
			}
			else if(qL[lIndex] == '-' || qL[lIndex] == '.'){
				//insert a gap into the right sequences
				qLNew += qL[lIndex];
				rLNew += rL[lIndex];
				lIndex++;
				
				if(rIndex != rLength){
					qRNew += '-';
					rRNew += '-';
				}
				else{
					qRNew += '.';
					rRNew += '.';
				}
			}
			else if(qR[rIndex] == '-' || qR[rIndex] == '.'){
				//insert a gap into the left sequences
				qRNew += qR[rIndex];
				rRNew += rR[rIndex];
				rIndex++;
				

				if(lIndex != lLength){
					qLNew += '-';
					rLNew += '-';
				}
				else{
					qLNew += '.';
					rLNew += '.';
				}
				
			}
		}
		
		qAlign = qLNew;
		aAlign = rLNew;
		bAlign = rRNew;
		
		bool qStart = 0;
		bool aStart = 0;
		bool bStart = 0;
		
		for(int i=0;i<qAlign.length();i++){
			
			if (m->getControl_pressed()) { return 0; }
			
			if(qStart == 0){
				if(qAlign[i] == '-')	{	qAlign[i] = '.';	}
				else					{	qStart = 1;			}
			}
			if(aStart == 0){
				if(aAlign[i] == '-')	{	aAlign[i] = '.';	}
				else					{	aStart = 1;			}
			}
			if(bStart == 0){
				if(bAlign[i] == '-')	{	bAlign[i] = '.';	}
				else					{	bStart = 1;			}
			}
			if(aStart == 1 && bStart == 1 && qStart == 1){
				break;
			}
		}
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "Perseus", "threeWayAlign");
		exit(1);
	}
}

/**************************************************************************************************/

double Perseus::calcLoonIndex(string query, string parent1, string parent2, int breakPoint, vector<vector<double> >& binMatrix){
	try {
		string queryAln, leftParentAln, rightParentAln;
		threeWayAlign(query, parent1, parent2, queryAln, leftParentAln, rightParentAln);
		
		int alignLength = queryAln.length();

		int endPos = alignLength;
		for(int i=alignLength-1;i>=0; i--){
			if(queryAln[i] != '.' && leftParentAln[i] != '.' && rightParentAln[i] != '.'){
				endPos = i + 1;
				break;
			}
		}
		
		int diffToLeftCount = 0;
		vector<int> diffToLeftMap(alignLength, 0);
		
		int diffToRightCount = 0;
		vector<int> diffToRightMap(alignLength, 0);
		
		for(int i=0;i<endPos;i++){
			
			if (m->getControl_pressed()) { return 0; }
			
			if(queryAln[i] != leftParentAln[i]){
				diffToLeftMap[diffToLeftCount] = i;
				diffToLeftCount++;
			}
			
			if(queryAln[i] != rightParentAln[i]){
				diffToRightMap[diffToRightCount] = i;
				diffToRightCount++;
			}
		}
		
		
		
		diffToLeftMap[diffToLeftCount] = endPos;
		diffToRightMap[diffToRightCount] = endPos;
		
		int indexL = 0;
		int indexR = 0;
		int indexS = 0;
		
		vector<int> diffs;
		vector<int> splits;
		
		splits.push_back(-1);
		diffs.push_back(diffToRightCount);
		indexS++;
			
		while(indexL < diffToLeftCount || indexR < diffToRightCount){
			
			if (m->getControl_pressed()) { return 0; }
			
			if(diffToLeftMap[indexL] <= diffToRightMap[indexR]){
				diffs.push_back(diffs[indexS - 1] + 1);
				splits.push_back(diffToLeftMap[indexL]);
				
				indexL++;
				indexS++;			
			}
			else if(diffToLeftMap[indexL] > diffToRightMap[indexR]) {
				diffs.push_back(diffs[indexS - 1] - 1);
				splits.push_back(diffToRightMap[indexR]);
				
				indexR++;
				indexS++;			
			}
		}
		
		int minDiff = PERSEUSMAXINT;
		int minIndex = -1;
		for(int i=0;i<indexS;i++){
			
			if (m->getControl_pressed()) { return 0; }
			
			if(diffs[i] < minDiff){
				minDiff = diffs[i];
				minIndex = i;
			}
		}
		
		int splitPos = endPos;
		if(minIndex < indexS - 1){
			splitPos = (splits[minIndex]+splits[minIndex+1]) / 2;
		}
		
		int diffToChimera = 0;
		int leftDiffToP1 = 0;
		int rightDiffToP1 = 0;
		int leftDiffToP2 = 0;
		int rightDiffToP2 = 0;
		
		for(int i=0;i<endPos;i++){
			
			if (m->getControl_pressed()) { return 0; }
			
			char bQuery = queryAln[i];
			char bP1 = leftParentAln[i];
			char bP2 = rightParentAln[i];
			
			char bConsensus = bQuery;
			if(bP1 == bP2){	bConsensus = bP1;	}
			
			if(bConsensus != bQuery){
				diffToChimera++;
			}
			
			if(bConsensus != bP1){
				if(i <= splitPos){
					leftDiffToP1++;
				}
				else{
					rightDiffToP1++;				
				}
			}
			if(bConsensus != bP2){
				if(i <= splitPos){
					leftDiffToP2++;
				}
				else{
					rightDiffToP2++;				
				}
			}
		}		
		

		int diffToClosestParent, diffToFurtherParent;
		int xA, xB, yA, yB;
		double aFraction, bFraction;
		
		if(diffToLeftCount <= diffToRightCount){	//if parent 1 is closer

			diffToClosestParent = leftDiffToP1 + rightDiffToP1;
			xA = leftDiffToP1;
			xB = rightDiffToP1;
			
			diffToFurtherParent = leftDiffToP2 + rightDiffToP2;
			yA = leftDiffToP2;
			yB = rightDiffToP2;
			
			aFraction = double(splitPos + 1)/(double) endPos;
			bFraction = 1 - aFraction;
			
		}
		else{												//if parent 2 is closer

			diffToClosestParent = leftDiffToP2 + rightDiffToP2;
			xA = rightDiffToP2;
			xB = leftDiffToP2;
			
			diffToFurtherParent = leftDiffToP1 + rightDiffToP1;
			yA = rightDiffToP1;
			yB = leftDiffToP1;
			
			bFraction = double(splitPos + 1)/(double) endPos;
			aFraction = 1 - bFraction;

		}
		
		double loonIndex = 0;
		
		int totalDifference = diffToClosestParent + diffToChimera;
		
		if(totalDifference > 0){
			double prob = 0;
			
			for(int i=diffToClosestParent;i<=totalDifference;i++){
				prob += binMatrix[totalDifference][i] * pow(0.50, i) * pow(0.50, totalDifference - i);
			}
			loonIndex += -log(prob);
		}
		
		if(diffToFurtherParent > 0){
			double prob = 0;
			
			for(int i=yA;i<=diffToFurtherParent;i++){
				prob += binMatrix[diffToFurtherParent][i] * pow(aFraction, i) * pow(1-aFraction, diffToFurtherParent - i);
			}
			loonIndex += -log(prob);
		}
		
		if(diffToClosestParent > 0){
			double prob = 0;
			
			for(int i=xB;i<=diffToClosestParent;i++){
				prob += binMatrix[diffToClosestParent][i] * pow(bFraction, i) * pow(1-bFraction, diffToClosestParent - i);
			}
			loonIndex += -log(prob);
		}
		
		return loonIndex;
	}
	catch(exception& e) {
		m->errorOut(e, "Perseus", "calcLoonIndex");
		exit(1);
	}
}

/**************************************************************************************************/

double Perseus::calcBestDistance(string query, string reference){
	try {
		int alignLength = query.length();
		int mismatch = 0;
		int counter = 0;
		
		for(int i=0;i<alignLength;i++){
			
			if (m->getControl_pressed()) { return 0; }
			
			if((query[i] != '.' || reference[i] != '.') && (query[i] != '-' && reference[i] != '-')){
				if(query[i] != reference[i]){	mismatch++;	}
				counter++;			
			}
		}
		
		return (double)mismatch / (double)counter;
	}
	catch(exception& e) {
		m->errorOut(e, "Perseus", "calcBestDistance");
		exit(1);
	}
}

/**************************************************************************************************/

double Perseus::classifyChimera(double singleDist, double cIndex, double loonIndex, double alpha, double beta){
	try {
		double difference = cIndex - singleDist;	//y
		double probability;
		
		if(cIndex >= 0.15 || difference > 0.00){
			probability = 0.0000;
		}
		else{
			probability = 1.0 / (1.0 + exp(-(alpha + beta * loonIndex)));
		}
		
		return probability;
	}
	catch(exception& e) {
		m->errorOut(e, "Perseus", "classifyChimera");
		exit(1);
	}
}
/**************************************************************************************************/
