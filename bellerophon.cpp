/*
 *  bellerophon.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 7/9/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "bellerophon.h"
#include "eachgapdist.h"
#include "ignoregaps.h"
#include "onegapdist.h"


//***************************************************************************************************************

Bellerophon::Bellerophon(string name, string o)  {
	try {
		fastafile = name;
		outputDir = o;
	}
	catch(exception& e) {
		errorOut(e, "Bellerophon", "Bellerophon");
		exit(1);
	}
}

//***************************************************************************************************************
void Bellerophon::print(ostream& out) {
	try {
		int above1 = 0;
		out << "Name\tScore\tLeft\tRight\t" << endl;
		//output prefenence structure to .chimeras file
		for (int i = 0; i < pref.size(); i++) {
			out << pref[i].name << '\t' << setprecision(3) << pref[i].score[0] << '\t' << pref[i].leftParent[0] << '\t' << pref[i].rightParent[0] << endl;
			
			//calc # of seqs with preference above 1.0
			if (pref[i].score[0] > 1.0) { 
				above1++; 
				mothurOut(pref[i].name + " is a suspected chimera at breakpoint " + toString(pref[i].midpoint)); mothurOutEndLine();
				mothurOut("It's score is " + toString(pref[i].score[0]) + " with suspected left parent " + pref[i].leftParent[0] + " and right parent " + pref[i].rightParent[0]); mothurOutEndLine();
			}
		}
		
		//output results to screen
		mothurOutEndLine();
		mothurOut("Sequence with preference score above 1.0: " + toString(above1)); mothurOutEndLine();
		int spot;
		spot = pref.size()-1;
		mothurOut("Minimum:\t" + toString(pref[spot].score[0])); mothurOutEndLine();
		spot = pref.size() * 0.975;
		mothurOut("2.5%-tile:\t" + toString(pref[spot].score[0])); mothurOutEndLine();
		spot = pref.size() * 0.75;
		mothurOut("25%-tile:\t" + toString(pref[spot].score[0])); mothurOutEndLine();
		spot = pref.size() * 0.50;
		mothurOut("Median: \t" + toString(pref[spot].score[0])); mothurOutEndLine();
		spot = pref.size() * 0.25;
		mothurOut("75%-tile:\t" + toString(pref[spot].score[0])); mothurOutEndLine();
		spot = pref.size() * 0.025;
		mothurOut("97.5%-tile:\t" + toString(pref[spot].score[0])); mothurOutEndLine();
		spot = 0;
		mothurOut("Maximum:\t" + toString(pref[spot].score[0])); mothurOutEndLine();

	}
	catch(exception& e) {
		errorOut(e, "Bellerophon", "print");
		exit(1);
	}
}

//********************************************************************************************************************
//sorts highest score to lowest
inline bool comparePref(Preference left, Preference right){
	return (left.score[0] > right.score[0]);	
}

//***************************************************************************************************************
int Bellerophon::getChimeras() {
	try {
		
		//do soft filter
		if (filter)  {
			string optionString = "fasta=" + fastafile + ", soft=50";
			if (outputDir != "") { optionString += ", outputdir=" + outputDir; }
			
			filterSeqs = new FilterSeqsCommand(optionString);
			filterSeqs->execute();
			delete filterSeqs;
			
			//reset fastafile to filtered file
			if (outputDir == "") { fastafile = getRootName(fastafile) + "filter.fasta"; }
			else				 { fastafile = outputDir + getRootName(getSimpleName(fastafile)) + "filter.fasta"; }
			
		}
		
		distCalculator = new eachGapDist();
		
		//read in sequences
		seqs = readSeqs(fastafile);
		
		if (unaligned) { mothurOut("Your sequences need to be aligned when you use the bellerophon method."); mothurOutEndLine(); return 1;  }
		
		int numSeqs = seqs.size();
		
		if (numSeqs == 0) { mothurOut("Error in reading you sequences."); mothurOutEndLine(); exit(1); }
		
		//set default window to 25% of sequence length
		string seq0 = seqs[0]->getAligned();
		if (window == 0) { window = seq0.length() / 4;  }
		else if (window > (seq0.length() / 2)) {  
			mothurOut("Your sequence length is = " + toString(seq0.length()) + ". You have selected a window size greater than the length of half your aligned sequence. I will run it with a window size of " + toString((seq0.length() / 2))); mothurOutEndLine();
			window = (seq0.length() / 2);
		}
		
		if (increment > (seqs[0]->getAlignLength() - (2*window))) { 
			if (increment != 10) {
			
				mothurOut("You have selected a increment that is too large. I will use the default."); mothurOutEndLine();
				increment = 10;
				if (increment > (seqs[0]->getAlignLength() - (2*window))) {  increment = 0;  }
				
			}else{ increment = 0; }
		}
		
		if (increment == 0) { iters = 1; }
		else { iters = ((seqs[0]->getAlignLength() - (2*window)) / increment); }
		
		//initialize pref
		pref.resize(numSeqs);  
		
		for (int i = 0; i < numSeqs; i++ ) { 
			pref[i].leftParent.resize(2); pref[i].rightParent.resize(2); pref[i].score.resize(2);   pref[i].closestLeft.resize(2); pref[i].closestRight.resize(3);
			pref[i].name = seqs[i]->getName();
			pref[i].score[0] = 0.0;  pref[i].score[1] = 0.0; 
			pref[i].closestLeft[0] = 100000.0;  pref[i].closestLeft[1] = 100000.0;  
			pref[i].closestRight[0] = 100000.0;  pref[i].closestRight[1] = 100000.0;  
		}

		int midpoint = window;
		int count = 0;
		while (count < iters) {
				
				//create 2 vectors of sequences, 1 for left side and one for right side
				vector<Sequence> left;  vector<Sequence> right;
				
				for (int i = 0; i < seqs.size(); i++) {
//cout << "midpoint = " << midpoint << "\twindow = " << window << endl;
//cout << "whole = " << seqs[i]->getAligned().length() << endl;
					//save left side
					string seqLeft = seqs[i]->getAligned().substr(midpoint-window, window);
					Sequence tempLeft;
					tempLeft.setName(seqs[i]->getName());
					tempLeft.setAligned(seqLeft);
					left.push_back(tempLeft);
//cout << "left = " << tempLeft.getAligned().length() << endl;			
					//save right side
					string seqRight = seqs[i]->getAligned().substr(midpoint, window);
					Sequence tempRight;
					tempRight.setName(seqs[i]->getName());
					tempRight.setAligned(seqRight);
					right.push_back(tempRight);
//cout << "right = " << seqRight.length() << endl;	
				}
				
				//adjust midpoint by increment
				midpoint += increment;
				
				
				//this should be parallelized
				//perference = sum of (| distance of my left to sequence j's left - distance of my right to sequence j's right | )
				//create a matrix containing the distance from left to left and right to right
				//calculate distances
				SparseMatrix* SparseLeft = new SparseMatrix();
				SparseMatrix* SparseRight = new SparseMatrix();
				
				createSparseMatrix(0, left.size(), SparseLeft, left);
				createSparseMatrix(0, right.size(), SparseRight, right);
				
				vector<SeqMap> distMapRight;
				vector<SeqMap> distMapLeft;
				
				// Create a data structure to quickly access the distance information.
				//this is from thallingers reimplementation on get.oturep
				// It consists of a vector of distance maps, where each map contains
				// all distances of a certain sequence. Vector and maps are accessed
				// via the index of a sequence in the distance matrix
				distMapRight = vector<SeqMap>(numSeqs); 
				distMapLeft = vector<SeqMap>(numSeqs); 
				//cout << "left" << endl << endl;
				for (MatData currentCell = SparseLeft->begin(); currentCell != SparseLeft->end(); currentCell++) {
					distMapLeft[currentCell->row][currentCell->column] = currentCell->dist;
					//cout << " i = " << currentCell->row << " j = " << currentCell->column << " dist = " << currentCell->dist << endl;
				}
				//cout << "right" << endl << endl;
				for (MatData currentCell = SparseRight->begin(); currentCell != SparseRight->end(); currentCell++) {
					distMapRight[currentCell->row][currentCell->column] = currentCell->dist;
					//cout << " i = " << currentCell->row << " j = " << currentCell->column << " dist = " << currentCell->dist << endl;
				}
				
				delete SparseLeft;
				delete SparseRight;
				
				//fill preference structure
				generatePreferences(distMapLeft, distMapRight, midpoint);
				
				count++;
				
		}
		
		delete distCalculator;
		
		//rank preference score to eachother
		float dme = 0.0;
		float expectedPercent = 1 / (float) (pref.size());
		
		for (int i = 0; i < pref.size(); i++) {	 dme += pref[i].score[0];  }
	
		for (int i = 0; i < pref.size(); i++) {

			//gives the actual percentage of the dme this seq adds
			pref[i].score[0] = pref[i].score[0] / dme;
			
			//how much higher or lower is this than expected
			pref[i].score[0] = pref[i].score[0] / expectedPercent;
		
		}
		
		//sort Preferences highest to lowest
		sort(pref.begin(), pref.end(), comparePref);
		
		return 0;

	}
	catch(exception& e) {
		errorOut(e, "Bellerophon", "getChimeras");
		exit(1);
	}
}

/***************************************************************************************************************/
int Bellerophon::createSparseMatrix(int startSeq, int endSeq, SparseMatrix* sparse, vector<Sequence> s){
	try {

		for(int i=startSeq; i<endSeq; i++){
			
			for(int j=0;j<i;j++){
			
				distCalculator->calcDist(s[i], s[j]);
				float dist = distCalculator->getDist();
			
				PCell temp(i, j, dist);
				sparse->addCell(temp);
				
			}
		}
		
		return 1;
	}
	catch(exception& e) {
		errorOut(e, "Bellerophon", "createSparseMatrix");
		exit(1);
	}
}
/***************************************************************************************************************/
void Bellerophon::generatePreferences(vector<SeqMap> left, vector<SeqMap> right, int mid){
	try {
		
		float dme = 0.0;
		SeqMap::iterator itR;
		SeqMap::iterator itL;
		
		//initialize pref[i]
		for (int i = 0; i < pref.size(); i++) {
			pref[i].score[1] = 0.0;
			pref[i].closestLeft[1] = 100000.0; 
			pref[i].closestRight[1] = 100000.0; 
			pref[i].leftParent[1] = "";
			pref[i].rightParent[1] = "";
		}
	
		for (int i = 0; i < left.size(); i++) {
			
			SeqMap currentLeft = left[i];    //example i = 3;   currentLeft is a map of 0 to the distance of sequence 3 to sequence 0,
												//										1 to the distance of sequence 3 to sequence 1,
												//										2 to the distance of sequence 3 to sequence 2.
			SeqMap currentRight = right[i];		// same as left but with distances on the right side.
			
			for (int j = 0; j < i; j++) {
				
				itL = currentLeft.find(j);
				itR = currentRight.find(j);
//cout << " i = " << i << " j = " << j << " distLeft = " << itL->second << endl;
//cout << " i = " << i << " j = " << j << " distright = " << itR->second << endl;
				
				//if you can find this entry update the preferences
				if ((itL != currentLeft.end()) && (itR != currentRight.end())) {
				
					if (!correction) {
						pref[i].score[1] += abs((itL->second - itR->second));
						pref[j].score[1] += abs((itL->second - itR->second));
//cout << "left " << i << " " << j << " = " << itL->second << " right " << i << " " << j << " = " << itR->second << endl;
//cout << "abs = " << abs((itL->second - itR->second)) << endl;
//cout << i << " score = " << pref[i].score[1] << endl;
//cout << j << " score = " << pref[j].score[1] << endl;
					}else {
						pref[i].score[1] += abs((sqrt(itL->second) - sqrt(itR->second)));
						pref[j].score[1] += abs((sqrt(itL->second) - sqrt(itR->second)));
//cout << "left " << i << " " << j << " = " << itL->second << " right " << i << " " << j << " = " << itR->second << endl;
//cout << "abs = " << abs((sqrt(itL->second) - sqrt(itR->second))) << endl;
//cout << i << " score = " << pref[i].score[1] << endl;
//cout << j << " score = " << pref[j].score[1] << endl;
					}
//cout << "pref[" << i << "].closestLeft[1] = "	<< 	pref[i].closestLeft[1] << " parent = " << pref[i].leftParent[1] << endl;			
					//are you the closest left sequence
					if (itL->second < pref[i].closestLeft[1]) {  

						pref[i].closestLeft[1] = itL->second;
						pref[i].leftParent[1] = seqs[j]->getName();
//cout << "updating closest left to " << pref[i].leftParent[1] << endl;
					}
//cout << "pref[" << j << "].closestLeft[1] = "	<< 	pref[j].closestLeft[1] << " parent = " << pref[j].leftParent[1] << endl;	
					if (itL->second < pref[j].closestLeft[1]) { 
						pref[j].closestLeft[1] = itL->second;
						pref[j].leftParent[1] = seqs[i]->getName();
//cout << "updating closest left to " << pref[j].leftParent[1] << endl;
					}
					
					//are you the closest right sequence
					if (itR->second < pref[i].closestRight[1]) {   
						pref[i].closestRight[1] = itR->second;
						pref[i].rightParent[1] = seqs[j]->getName();
					}
					if (itR->second < pref[j].closestRight[1]) {   
						pref[j].closestRight[1] = itR->second;
						pref[j].rightParent[1] = seqs[i]->getName();
					}
					
				}
			}
		
		}
		
		
		  
		//calculate the dme
		int count0 = 0;
		for (int i = 0; i < pref.size(); i++) {	 dme += pref[i].score[1];  if (pref[i].score[1] == 0.0) { count0++; }  }
		
		float expectedPercent = 1 / (float) (pref.size() - count0);
//cout << endl << "dme = " << dme << endl;
		//recalculate prefernences based on dme
		for (int i = 0; i < pref.size(); i++) {
//cout << "unadjusted pref " << i << " = " << pref[i].score[1] << endl;	
			// gives the actual percentage of the dme this seq adds
			pref[i].score[1] = pref[i].score[1] / dme;
			
			//how much higher or lower is this than expected
			pref[i].score[1] = pref[i].score[1] / expectedPercent;
			
			//pref[i].score[1] = dme / (dme - 2 * pref[i].score[1]);
			
			//so a non chimeric sequence would be around 1, and a chimeric would be signifigantly higher.
//cout << "adjusted pref " << i << " = " << pref[i].score[1] << endl;					
		}
		
		//is this score bigger then the last score
		for (int i = 0; i < pref.size(); i++) {	 
			
			//update biggest score
			if (pref[i].score[1] > pref[i].score[0]) {
				pref[i].score[0] = pref[i].score[1];
				pref[i].leftParent[0] = pref[i].leftParent[1];
				pref[i].rightParent[0] = pref[i].rightParent[1];
				pref[i].closestLeft[0] = pref[i].closestLeft[1];
				pref[i].closestRight[0] = pref[i].closestRight[1];
				pref[i].midpoint = mid;
			}
			
		}

	}
	catch(exception& e) {
		errorOut(e, "Bellerophon", "generatePreferences");
		exit(1);
	}
}
/**************************************************************************************************/

