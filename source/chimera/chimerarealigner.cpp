/*
 *  chimerarealigner.cpp
 *  Mothur
 *
 *  Created by westcott on 2/12/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "chimerarealigner.h"
#include "needlemanoverlap.hpp"
#include "nast.hpp"

//***************************************************************************************************************
ChimeraReAligner::ChimeraReAligner()  {  m = MothurOut::getInstance(); }
//***************************************************************************************************************
ChimeraReAligner::~ChimeraReAligner() = default;	
//***************************************************************************************************************
void ChimeraReAligner::reAlign(Sequence* query, vector<string> parents) {
	try {

		if(parents.size() != 0){

			alignmentLength = query->getAlignLength();	//x
			int queryUnalignedLength = query->getNumBases();	//y
			
			buildTemplateProfile(parents);
			
			createAlignMatrix(queryUnalignedLength, alignmentLength);
			fillAlignMatrix(query->getUnaligned());
			query->setAligned(getNewAlignment(query->getUnaligned()));
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraReAligner", "reAlign");
		exit(1);
	}
}
/***************************************************************************************************************/

void ChimeraReAligner::buildTemplateProfile(vector<string> parents) {
	try{	
		int numParents = parents.size();

		profile.resize(alignmentLength);
				
		for(int i=0;i<numParents;i++){
			string seq = parents[i];

			for(int j=0;j<alignmentLength;j++){
					
				
				if(seq[j] == 'A')		{	profile[j].A++;		}
				else if(seq[j] == 'T')	{	profile[j].T++;		}
				else if(seq[j] == 'G')	{	profile[j].G++;		}
				else if(seq[j] == 'C')	{	profile[j].C++;		}
				else if(seq[j] == '-')	{	profile[j].Gap++;	}
				else if(seq[j] == '.')	{	profile[j].Gap++;	}
				
			}
		}
		

		for(int i=0;i<alignmentLength;i++){
			profile[i].Chars = profile[i].A + profile[i].T + profile[i].G + profile[i].C;
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraReAligner", "buildTemplateProfile");
		exit(1);
	}
}
	
 
/***************************************************************************************************************/

void ChimeraReAligner::createAlignMatrix(int queryUnalignedLength, int alignmentLength){
	
	try{
		alignMatrix.resize(alignmentLength+1);
		for(int i=0;i<=alignmentLength;i++){
			alignMatrix[i].resize(queryUnalignedLength+1);
		}

		for(int i=1;i<=alignmentLength;i++)		{	alignMatrix[i][0].direction = 'l';	}
		for(int j=1;j<=queryUnalignedLength;j++){	alignMatrix[0][j].direction = 'u';	}
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraReAligner", "createAlignMatrix");
		exit(1);
	}
}

/***************************************************************************************************************/

void ChimeraReAligner::fillAlignMatrix(string query){
	try{
		int GAP = -4;
		
		int nrows = alignMatrix.size()-1;
		int ncols = alignMatrix[0].size()-1;
				
		for(int i=1;i<=nrows;i++){
			
			bases p = profile[i-1];
			int numChars = p.Chars;
			
			for(int j=1;j<=ncols;j++){
			
				char q = query[j-1];
				
				//	score it for if there was a match
				int maxScore = calcMatchScore(p, q) + alignMatrix[i-1][j-1].score;
				int maxDirection = 'd';
				
				//	score it for if there was a gap in the query
				int score = alignMatrix[i-1][j].score + (numChars * GAP);
				if (score > maxScore) {
					maxScore = score;
					maxDirection = 'l';
				}
				
				alignMatrix[i][j].score = maxScore;
				alignMatrix[i][j].direction = maxDirection;
				
			}
		}
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraReAligner", "fillAlignMatrix");
		exit(1);
	}
}

/***************************************************************************************************************/

int ChimeraReAligner::calcMatchScore(bases p, char q){
	try{
		
		int MATCH = 5;
		int MISMATCH = -4;
		
		int score = 0;
		
		if(q == 'G')		{	score = (MATCH * p.G + MISMATCH * (p.A + p.T + p.C + p.Gap));		}
		else if(q == 'A')	{	score = (MATCH * p.A + MISMATCH * (p.G + p.T + p.C + p.Gap));		}
		else if(q == 'T')	{	score = (MATCH * p.T + MISMATCH * (p.G + p.A + p.C + p.Gap));		}
		else if(q == 'C')	{	score = (MATCH * p.C + MISMATCH * (p.G + p.A + p.T + p.Gap));		}
		else				{	score = (MATCH * p.A + MISMATCH * (p.G + p.T + p.C + p.Gap));		}

		return score;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraReAligner", "calcMatchScore");
		exit(1);
	}
}

/***************************************************************************************************************/

string ChimeraReAligner::getNewAlignment(string query){
	try{
		string queryAlignment(alignmentLength, '.');
		string referenceAlignment(alignmentLength, '.');
		
		
		int maxScore = -99999999;
		
		int nrows = alignMatrix.size()-1;
		int ncols = alignMatrix[0].size()-1;

		int bestCol = -1;
		int bestRow = -1;
		
		for(int i=1;i<=nrows;i++){
			int score = alignMatrix[i][ncols].score;
			if (score > maxScore) {
				maxScore = score;
				bestRow = i;
				bestCol = ncols;
			}
		}
		
		for(int j=1;j<=ncols;j++){
			int score = alignMatrix[nrows][j].score;
			if (score > maxScore) {
				maxScore = score;
				bestRow = nrows;
				bestCol = j;
			}
		}
		
		int currentRow = bestRow;
		int currentCol = bestCol;
		
		int alignmentPosition = 0;
		if(currentRow < alignmentLength){
			for(int i=alignmentLength;i>currentRow;i--){
				alignmentPosition++;
			}
		}
		
		AlignCell c = alignMatrix[currentRow][currentCol];
		while(c.direction != 'x'){
			
			char q;

			if(c.direction == 'd'){
				q = query[currentCol-1];
				currentCol--;
				currentRow--;
			}
			
			
			else if (c.direction == 'u') {
				break;
			}					
			else if(c.direction == 'l'){
				char gapChar;
				if(currentCol == 0)	{	gapChar = '.';	}
				else				{	gapChar = '-';	}
								
				q = gapChar;
				currentRow--;
			}
            else{ m->mothurOut("[ERROR]: Unexpected case in ChimeraReAligner::getNewAlignment, aborting.\n"); m->setControl_pressed(true); }

			queryAlignment[alignmentPosition] = q;
			alignmentPosition++;
			c = alignMatrix[currentRow][currentCol];
		}

//		need to reverse the string
		string flipSeq = "";
		for(int i=alignmentLength-1;i>=0;i--){
			flipSeq += queryAlignment[i];			
		}

		return flipSeq;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraReAligner", "getNewAlignment");
		exit(1);
	}
}

/***************************************************************************************************************/

// Sequence* ChimeraReAligner::getSequence(string name) {
//	try{
//		Sequence* temp;
//		
//		//look through templateSeqs til you find it
//		int spot = -1;
//		for (int i = 0; i < templateSeqs.size(); i++) {
//			if (name == templateSeqs[i]->getName()) {  
//				spot = i;
//				break;
//			}
//		}
//		
//		if(spot == -1) { m->mothurOut("Error: Could not find sequence.\n");  return nullptr; }
//		
//		temp = new Sequence(templateSeqs[spot]->getName(), templateSeqs[spot]->getAligned());
//		
//		return temp;
//	}
//	catch(exception& e) {
//		m->errorOut(e, "ChimeraReAligner", "getSequence");
//		exit(1);
//	}
//}

//***************************************************************************************************************/
