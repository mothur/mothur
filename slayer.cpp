/*
 *  slayer.cpp
 *  Mothur
 *
 *  Created by westcott on 9/25/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */

#include "slayer.h"

/***********************************************************************/
Slayer::Slayer(int win, int increment, int parentThreshold, float div) :
		windowSize(win), windowStep(increment), parentFragmentThreshold(parentThreshold), divRThreshold(div) {}
/***********************************************************************/
void Slayer::getResults(Sequence* query, vector<Sequence*> refSeqs) {
	try {
		
		for (int i = 0; i < refSeqs.size(); i++) {
		
			for (int j = i+1; j < refSeqs.size(); j++) {
			
				//make copies of query and each parent because runBellerophon removes gaps and messes them up
				Sequence* q = new Sequence(query->getName(), query->getAligned());
				Sequence* leftParent = new Sequence(refSeqs[i]->getName(), refSeqs[i]->getAligned());
				Sequence* rightParent = new Sequence(refSeqs[j]->getName(), refSeqs[j]->getAligned());
				
				vector<data_struct> divs = runBellerophon(q, leftParent, rightParent);
				
				delete q;
				delete leftParent;
				delete rightParent;
			}
			
		}
		
		
	}
	catch(exception& e) {
		errorOut(e, "Slayer", "getResults");
		exit(1);
	}
}
/***********************************************************************/
vector<data_struct> Slayer::runBellerophon(Sequence* query, Sequence* parentA, Sequence* parentB) {
	try{
		
		vector<data_struct> data;
		
		//vertical filter
		vector<Sequence*> temp;
		verticalFilter(temp);
		
		int alignLength = query->getAligned().length();
		
		
		
		
		return data;
		
	}
	catch(exception& e) {
		errorOut(e, "Slayer", "runBellerophon");
		exit(1);
	}
}
/***********************************************************************/
float Slayer::computePercentID(string queryFrag, string parent, int left, int right) {
	try {
		int total = 0;
		int matches = 0;
	
		for (int i = left; i <= right; i++) {
			total++;
			if (queryFrag[i] == parent[i]) {
				matches++;
			}
		}

		float percentID =( matches/(float)total) * 100;
		
		return percentID;
	}
	catch(exception& e) {
		errorOut(e, "Slayer", "computePercentID");
		exit(1);
	}
}
/***********************************************************************/
//this is a vertical filter
void Slayer::verticalFilter(vector<Sequence*> seqs) {
	try {
		vector<int> gaps;	gaps.resize(seqs[0]->getAligned().length(), 0);
		
		string filterString = (string(seqs[0]->getAligned().length(), '1'));
		
		//for each sequence
		for (int i = 0; i < seqs.size(); i++) {
		
			string seqAligned = seqs[i]->getAligned();
			
			for (int j = 0; j < seqAligned.length(); j++) {
				//if this spot is a gap
				if ((seqAligned[j] == '-') || (seqAligned[j] == '.'))	{	gaps[j]++;	}
			}
		}
		
		//zero out spot where all sequences have blanks
		int numColRemoved = 0;
		for(int i = 0; i < seqs[0]->getAligned().length(); i++){
			if(gaps[i] == seqs.size())	{	filterString[i] = '0'; 	numColRemoved++;  }
		}
		
		//for each sequence
		for (int i = 0; i < seqs.size(); i++) {
		
			string seqAligned = seqs[i]->getAligned();
			string newAligned = "";
			
			for (int j = 0; j < seqAligned.length(); j++) {
				//if this spot is not a gap
				if (filterString[j] == '1') { newAligned += seqAligned[j]; }
			}
			
			seqs[i]->setAligned(newAligned);
		}
	}
	catch(exception& e) {
		errorOut(e, "Slayer", "verticalFilter");
		exit(1);
	}
}
/***********************************************************************/
