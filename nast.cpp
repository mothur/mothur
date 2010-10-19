/*
 *  nast.cpp
 *  
 *
 *  Created by Pat Schloss on 12/17/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 *	This is my implementation of the NAST (nearest alignment space termination) algorithm as described in:
 *
 *	DeSantis TZ, Hugenholtz P, Keller K, Brodie EL, Larsen N, Piceno YM, Phan R, & Anderson GL.  2006.  NAST: a multiple
 *		sequence alignment server for comparative analysis of 16S rRNA genes.  Nucleic Acids Research.  34:W394-9.
 *
 *	To construct an object one needs to provide a method of getting a pairwise alignment (alignment) and the template
 *	and candidate sequence that are to be aligned to each other.
 *
 */

#include "sequence.hpp"
#include "alignment.hpp"
#include "nast.hpp"

/**************************************************************************************************/

Nast::Nast(Alignment* method, Sequence* cand, Sequence* temp) : alignment(method), candidateSeq(cand), templateSeq(temp) {
	try {
		m = MothurOut::getInstance();
		maxInsertLength = 0;
	
		pairwiseAlignSeqs();	//	This is part A in Fig. 2 of DeSantis et al.
		regapSequences();		//	This is parts B-F in Fig. 2 of DeSantis et al.
		
	}
	catch(exception& e) {
		m->errorOut(e, "Nast", "Nast");
		exit(1);
	}
}

/**************************************************************************************************/

void Nast::pairwiseAlignSeqs(){	//	Here we call one of the pairwise alignment methods to align our unaligned candidate
								//	and template sequences
	try {	
		alignment->align(candidateSeq->getUnaligned(), templateSeq->getUnaligned());

		string candAln = alignment->getSeqAAln();
		string tempAln = alignment->getSeqBAln();

		if(candAln == ""){

			candidateSeq->setPairwise("");
			templateSeq->setPairwise(templateSeq->getUnaligned());

		}
		else{
			if(tempAln[0] == '-'){
				int pairwiseAlignmentLength = tempAln.length();	//	we need to make sure that the candidate sequence alignment
				for(int i=0;i<pairwiseAlignmentLength;i++){		//	starts where the template sequence alignment starts, if it
					if(isalpha(tempAln[i])){					//	starts early, we nuke the beginning of the candidate
						candAln = candAln.substr(i);			//	sequence
						tempAln = tempAln.substr(i);
						break;
					}
				}
			}
			int pairwiseAlignmentLength = tempAln.length();
			if(tempAln[pairwiseAlignmentLength-1] == '-'){		//	we need to make sure that the candidate sequence alignment
				for(int i=pairwiseAlignmentLength-1; i>=0; i--){//	ends where the template sequence alignment ends, if it runs
					if(isalpha(tempAln[i])){					//	long, we nuke the end of the candidate sequence
						candAln = candAln.substr(0,i+1);
						tempAln = tempAln.substr(0,i+1);
						break;
					}		
				}
			}

		}

		candidateSeq->setPairwise(candAln);					//	set the pairwise sequences in the Sequence objects for
		templateSeq->setPairwise(tempAln);					//	the candidate and template sequences
	}
	catch(exception& e) {
		m->errorOut(e, "Nast", "pairwiseAlignSeqs");
		exit(1);
	}	
}

/**************************************************************************************************/

void Nast::removeExtraGaps(string& candAln, string tempAln, string newTemplateAlign){

//	here we do steps C-F of Fig. 2 from DeSantis et al.
	try {
	
		//cout << candAln << endl;
		//cout << tempAln << endl;
		//cout << newTemplateAlign << endl;
		//cout << endl;
		
		int longAlignmentLength = newTemplateAlign.length();	
	
		for(int i=0; i<longAlignmentLength; i++){				//	use the long alignment as the standard
			int rightIndex, rightRoom, leftIndex, leftRoom;
	
			//	Part C of Fig. 2 from DeSantis et al.
			if((isalpha(newTemplateAlign[i]) != isalpha(tempAln[i]))){	//if there is a discrepancy between the regapped
				
				rightRoom = 0; leftRoom = 0;
				
				//	Part D of Fig. 2 from DeSantis et al.		//	template sequence and the official template sequence
				for(leftIndex=i-1;leftIndex>0;leftIndex--){		//	then we've got problems...
					if(!isalpha(candAln[leftIndex])){
						leftRoom = 1;	//count how far it is to the nearest gap on the LEFT side of the anomaly
						while(leftIndex-leftRoom>=0 && !isalpha(candAln[leftIndex-leftRoom]))	{	leftRoom++;		}
						break;
					}
				}

				for(rightIndex=i+1;rightIndex<longAlignmentLength-1;rightIndex++){
					if(!isalpha(candAln[rightIndex])){
						rightRoom = 1;	//count how far it is to the nearest gap on the RIGHT side of the anomaly
						while(rightIndex+rightRoom<longAlignmentLength && !isalpha(candAln[rightIndex+rightRoom]))	{	rightRoom++;	}
						break;
					}
				}
								
				int insertLength = 0;							//	figure out how long the anomaly is
				while(!isalpha(newTemplateAlign[i + insertLength]))	{	insertLength++;	}
				if(insertLength > maxInsertLength){	maxInsertLength = insertLength;	}
		
				if((leftRoom + rightRoom) >= insertLength){
	
					//	Parts D & E from Fig. 2 of DeSantis et al.
					if((i-leftIndex) <= (rightIndex-i)){		//	the left gap is closer - > move stuff left there's
	
						if(leftRoom >= insertLength){			//	enough room to the left to move
			//cout << "lr newTemplateAlign = " << newTemplateAlign.length() << '\t' << i << '\t' << insertLength << endl;
							string leftTemplateString = newTemplateAlign.substr(0,i);
							string rightTemplateString = newTemplateAlign.substr((i+insertLength));
							newTemplateAlign = leftTemplateString + rightTemplateString;
							longAlignmentLength = newTemplateAlign.length();
			//cout << "lr candAln = " << candAln.length() << '\t' << leftIndex << '\t'  << endl;				
							string leftCandidateString = candAln.substr(0,(leftIndex-insertLength+1));
							string rightCandidateString = candAln.substr((leftIndex+1));
							candAln = leftCandidateString + rightCandidateString;
		
						}
						else{									//	not enough room to the left, have to steal some space to
						
			//cout << "in else lr newTemplateAlign = " << newTemplateAlign.length() << '\t' << i << '\t' << insertLength << endl;
							string leftTemplateString = newTemplateAlign.substr(0,i);	//	the right
							string rightTemplateString = newTemplateAlign.substr((i+insertLength));
							newTemplateAlign = leftTemplateString + rightTemplateString;
							longAlignmentLength = newTemplateAlign.length();
			//cout << " in else lr candAln = " << candAln.length() << '\t' << " leftIndex = " << leftIndex << " leftroom = " << leftRoom << " rightIndex = " << rightIndex << '\t' << endl;					
							string leftCandidateString = candAln.substr(0,(leftIndex-leftRoom+1));
							string insertString = candAln.substr((leftIndex+1),(rightIndex-leftIndex-1));
							string rightCandidateString = candAln.substr((rightIndex+(insertLength-leftRoom)));
							candAln = leftCandidateString + insertString + rightCandidateString;
				
						}
					}
					else{										//	the right gap is closer - > move stuff right there's
						if(rightRoom >= insertLength){			//	enough room to the right to move
			//cout << "rr newTemplateAlign = " << newTemplateAlign.length() << '\t' << i << '\t' << i+insertLength << endl;
							string leftTemplateString = newTemplateAlign.substr(0,i);
							string rightTemplateString = newTemplateAlign.substr((i+insertLength));
							newTemplateAlign = leftTemplateString + rightTemplateString;
							longAlignmentLength = newTemplateAlign.length();
			//cout << "rr candAln = " << candAln.length() << '\t' << i << '\t' << rightIndex << '\t' << rightIndex+insertLength << endl;				
							string leftCandidateString = candAln.substr(0,rightIndex);
							string rightCandidateString = candAln.substr((rightIndex+insertLength));
							candAln = leftCandidateString + rightCandidateString;	
									
						}
						else{									//	not enough room to the right, have to steal some 	
							//	space to the left lets move left and then right...
					//cout << "in else rr newTemplateAlign = " << newTemplateAlign.length() << '\t' << i << '\t' << i+insertLength << endl;
							string leftTemplateString = newTemplateAlign.substr(0,i);
							string rightTemplateString = newTemplateAlign.substr((i+insertLength));
							newTemplateAlign = leftTemplateString + rightTemplateString;
							longAlignmentLength = newTemplateAlign.length();
					//cout << "in else rr candAln = " << candAln.length() << '\t' << '\t' << (leftIndex-(insertLength-rightRoom)+1) << '\t' <<  (leftIndex+1,rightIndex-leftIndex-1) << '\t' << (rightIndex+rightRoom) << endl;				
							string leftCandidateString = candAln.substr(0,(leftIndex-(insertLength-rightRoom)+1));
							string insertString = candAln.substr((leftIndex+1),(rightIndex-leftIndex-1));
							string rightCandidateString = candAln.substr((rightIndex+rightRoom));
							candAln = leftCandidateString + insertString + rightCandidateString;	
									
						}
					}
					i -= insertLength;

				}
				else{
			//	there could be a case where there isn't enough room in either direction to move stuff
//cout << "in else else newTemplateAlign = " << newTemplateAlign.length() << '\t' << i << '\t' << (i+leftRoom+rightRoom) << endl;
					string leftTemplateString = newTemplateAlign.substr(0,i);	
					string rightTemplateString = newTemplateAlign.substr((i+leftRoom+rightRoom));
					newTemplateAlign = leftTemplateString + rightTemplateString;
					longAlignmentLength = newTemplateAlign.length();
							
		//cout << "in else else newTemplateAlign = " << candAln.length() << '\t' << (leftIndex-leftRoom+1) << '\t' << (leftIndex+1) << '\t' << (rightIndex-leftIndex-1) << '\t' << (rightIndex+rightRoom) << endl;	
					string leftCandidateString = candAln.substr(0,(leftIndex-leftRoom+1));
					string insertString = candAln.substr((leftIndex+1),(rightIndex-leftIndex-1));
					string rightCandidateString = candAln.substr((rightIndex+rightRoom));
					candAln = leftCandidateString + insertString + rightCandidateString;
					
					i -= (leftRoom + rightRoom);
				}
			
//				i -= insertLength;
				
				//if i is negative, we want to remove the extra gaps to the right
				if (i < 0) { cout << "i is negative" << endl; }
			} 
		}
	}
	catch(exception& e) {
		m->errorOut(e, "Nast", "removeExtraGaps");
		exit(1);
	}	
}

/**************************************************************************************************/

void Nast::regapSequences(){	//This is essentially part B in Fig 2. of DeSantis et al.
	try { 
	//cout << candidateSeq->getName() << endl;
		string candPair = candidateSeq->getPairwise();
		string candAln = "";
		
		string tempPair = templateSeq->getPairwise();
		string tempAln = templateSeq->getAligned();		//	we use the template aligned sequence as our guide
		
		int pairwiseLength = candPair.length();
		int fullAlignLength = tempAln.length();
		
		if(candPair == ""){
			for(int i=0;i<fullAlignLength;i++)	{	candAln += '.';		}
			candidateSeq->setAligned(candAln);
			return;
		}
	
		int fullAlignIndex = 0;
		int pairwiseAlignIndex = 0;
		string newTemplateAlign = "";					//	this is going to be messy so we want a temporary template
		//	alignment string
		while(tempAln[fullAlignIndex] == '.' || tempAln[fullAlignIndex]  == '-'){
			candAln += '.';								//	add the initial '-' and '.' to the candidate and template
			newTemplateAlign += tempAln[fullAlignIndex];//	pairwise sequences
			fullAlignIndex++;
		}

		string lastLoop = "";
		
		while(pairwiseAlignIndex<pairwiseLength){
	//cout << pairwiseAlignIndex << '\t' << fullAlignIndex << '\t' << pairwiseLength << endl;
			if(isalpha(tempPair[pairwiseAlignIndex]) && isalpha(tempAln[fullAlignIndex])
			   && isalpha(candPair[pairwiseAlignIndex])){
				//  the template and candidate pairwise and template aligned have characters
				//	need to add character onto the candidatSeq.aligned sequence
				
				candAln += candPair[pairwiseAlignIndex];
				newTemplateAlign += tempPair[pairwiseAlignIndex];//
				
				pairwiseAlignIndex++;
				fullAlignIndex++;
			}
			else if(isalpha(tempPair[pairwiseAlignIndex]) && !isalpha(tempAln[fullAlignIndex])
					&& isalpha(candPair[pairwiseAlignIndex])){
				//	the template pairwise and candidate pairwise are characters and the template aligned is a gap
				//	need to insert gaps into the candidateSeq.aligned sequence
				
				candAln += '-';
				newTemplateAlign += '-';//
				fullAlignIndex++;
			}
			else if(!isalpha(tempPair[pairwiseAlignIndex]) && isalpha(tempAln[fullAlignIndex])
					&& isalpha(candPair[pairwiseAlignIndex])){
				//  the template pairwise is a gap and the template aligned and pairwise sequences have characters
				//	this is the alpha scenario.  add character to the candidateSeq.aligned sequence without progressing
				//	further through the tempAln sequence.
				
				candAln += candPair[pairwiseAlignIndex];
				newTemplateAlign += '-';//
				pairwiseAlignIndex++;
			}
			else if(isalpha(tempPair[pairwiseAlignIndex]) && isalpha(tempAln[fullAlignIndex])
					&& !isalpha(candPair[pairwiseAlignIndex])){
				//  the template pairwise and full alignment are characters and the candidate sequence has a gap
				//	should not be a big deal, just add the gap position to the candidateSeq.aligned sequence;
				
				candAln += candPair[pairwiseAlignIndex];
				newTemplateAlign += tempAln[fullAlignIndex];//
				fullAlignIndex++;			
				pairwiseAlignIndex++;
			}
			else if(!isalpha(tempPair[pairwiseAlignIndex]) && !isalpha(tempAln[fullAlignIndex])
					&& isalpha(candPair[pairwiseAlignIndex])){
				//	the template pairwise and aligned are gaps while the candidate pairwise has a character
				//	this would be an insertion, go ahead and add the character->seems to be the opposite of the alpha scenario
				
				candAln += candPair[pairwiseAlignIndex];
				newTemplateAlign += tempAln[fullAlignIndex];//
				pairwiseAlignIndex++;
				fullAlignIndex++;			
			}
			else if(isalpha(tempPair[pairwiseAlignIndex]) && !isalpha(tempAln[fullAlignIndex])
					&& !isalpha(candPair[pairwiseAlignIndex])){
				//	template pairwise has a character, but its full aligned sequence and candidate sequence have gaps
				//	this would happen like we need to add a gap.  basically the opposite of the alpha situation
				
				newTemplateAlign += tempAln[fullAlignIndex];//
				candAln += "-";
				fullAlignIndex++;			
			}
			else if(!isalpha(tempPair[pairwiseAlignIndex]) && isalpha(tempAln[fullAlignIndex])
					&& !isalpha(candPair[pairwiseAlignIndex])){
				//	template and candidate pairwise are gaps and the template aligned is not a gap this should not be possible
				//	would skip the gaps and not progress through full alignment sequence
				//	not tested yet
				
				m->mothurOut("We're into D " + toString(fullAlignIndex) + " " +  toString(pairwiseAlignIndex)); m->mothurOutEndLine();
				pairwiseAlignIndex++;
			}
			else{
				//	everything has a gap - not possible
				//	not tested yet
				
				m->mothurOut("We're into F " +  toString(fullAlignIndex) + " " +  toString(pairwiseAlignIndex)); m->mothurOutEndLine();
				pairwiseAlignIndex++;
				fullAlignIndex++;			
			}		
		}
		
		for(int i=fullAlignIndex;i<fullAlignLength;i++){
			candAln += '.';
			newTemplateAlign += tempAln[i];//
		}
		
		int start = 0;
		int end = candAln.length()-1;

		for(int i=0;i<candAln.length();i++){
			if(candAln[i] == 'Z' || !isalnum(candAln[i]))	{	candAln[i] = '.';	}	//	if we padded the alignemnt from
			else{			start = i;			break;		}							//	blast with Z's, change them to
		}																				//	'.' characters
		
		for(int i=candAln.length()-1;i>=0;i--){											//	ditto.
			if(candAln[i] == 'Z' || !isalnum(candAln[i]))	{	candAln[i] = '.';	}
			else{			end = i;			break;		}
		}
		
		for(int i=start;i<=end;i++){					//	go through the candidate alignment sequence and make sure that
			candAln[i] = toupper(candAln[i]);			//	everything is upper case
		}
		

		if(candAln.length() != tempAln.length()){		//	if the regapped candidate sequence is longer than the official
			removeExtraGaps(candAln, tempAln, newTemplateAlign);//	template alignment then we need to do steps C-F in Fig.
		}												//	2 of Desantis et al.

		candidateSeq->setAligned(candAln);
	//cout << "here" << endl;
	}
	catch(exception& e) {
		m->errorOut(e, "Nast", "regapSequences");
		exit(1);
	}	
}

/**************************************************************************************************/

float Nast::getSimilarityScore(){
	try {
	
		string cand = candidateSeq->getAligned();
		string temp = templateSeq->getAligned();
		int alignmentLength = temp.length();
		int mismatch = 0;
		int denominator = 0;
		
		for(int i=0;i<alignmentLength;i++){
			if(cand[i] == '-' && temp[i] == '-'){
				
			}
			else if(cand[i] != '.' && temp[i] != '.'){
				denominator++;
				
				if(cand[i] != temp[i]){
					mismatch++;
				}
			}
		}
		float similarity = 100 * (1. - mismatch / (float)denominator);
		if(denominator == 0){	similarity = 0.0000;	}
		
		return similarity;
		
	}
	catch(exception& e) {
		m->errorOut(e, "Nast", "getSimilarityScore");
		exit(1);
	}	
}

/**************************************************************************************************/

int Nast::getMaxInsertLength(){
	
	return maxInsertLength;
	
}
	
/**************************************************************************************************/
