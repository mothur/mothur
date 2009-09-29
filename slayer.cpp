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
string Slayer::getResults(Sequence* query, vector<Sequence*> refSeqs) {
	try {
cout << "refSeqs = " << refSeqs.size() << endl;		
		vector<data_struct> all; all.clear();
		
		for (int i = 0; i < refSeqs.size(); i++) {
		
			for (int j = i+1; j < refSeqs.size(); j++) {
			
				//make copies of query and each parent because runBellerophon removes gaps and messes them up
				Sequence* q = new Sequence(query->getName(), query->getAligned());
				Sequence* leftParent = new Sequence(refSeqs[i]->getName(), refSeqs[i]->getAligned());
				Sequence* rightParent = new Sequence(refSeqs[j]->getName(), refSeqs[j]->getAligned());
				
				vector<data_struct> divs = runBellerophon(q, leftParent, rightParent);
				
				vector<data_struct> selectedDivs;
				for (int k = 0; k < divs.size(); k++) {
				
					vector<snps> snpsLeft = getSNPS(divs[k].parentA.getAligned(), divs[k].querySeq.getAligned(), divs[k].parentB.getAligned(), divs[k].winLStart, divs[k].winLEnd);
					vector<snps> snpsRight = getSNPS(divs[k].parentA.getAligned(), divs[k].querySeq.getAligned(), divs[k].parentB.getAligned(), divs[k].winRStart, divs[k].winREnd);
					
					int numSNPSLeft = snpsLeft.size();
					int numSNPSRight = snpsRight.size();
					
					//require at least 3 SNPs on each side of the break
					if ((numSNPSLeft >= 3) && (numSNPSRight >= 3)) {
						
						int winSizeLeft = divs[k].winLEnd - divs[k].winLStart + 1;
						int winSizeRight = divs[k].winREnd - divs[k].winRStart + 1;
						
						float snpRateLeft = numSNPSLeft / (float) winSizeLeft;
						float snpRateRight = numSNPSRight / (float) winSizeRight;
						
						float logR = log(snpRateLeft / snpRateRight) / log(2);
						
						// do not accept excess snp ratio on either side of the break
						if (abs(logR) < 1 ) {  
							
							float BS_A, BS_B;
							bootstrapSNPS(snpsLeft, snpsRight, BS_A, BS_B);
						
							divs[k].bsa = BS_A;
							divs[k].bsb = BS_B;
						
							divs[k].bsMax = max(BS_A, BS_B);
						
							divs[k].chimeraMax = max(divs[k].qla_qrb, divs[k].qlb_qra);
						
							selectedDivs.push_back(divs[k]);
						}
					}
				}
				
				//save selected
				for (int m = 0; m < selectedDivs.size(); m++) {  all.push_back(selectedDivs[m]);	}
				
				delete q;
				delete leftParent;
				delete rightParent;
			}
		}
		
		
		// compute bootstrap support
		if (all.size() > 0) {
			//sort them
			sort(all.begin(), all.end(), compareDataStruct);
			reverse(all.begin(), all.end());
			
			outputResults = all;
			return "yes"; 
		}
		else {
			outputResults = all;
			return "no";
		}
	}
	catch(exception& e) {
		errorOut(e, "Slayer", "getResults");
		exit(1);
	}
}
/***********************************************************************/
vector<data_struct> Slayer::runBellerophon(Sequence* q, Sequence* pA, Sequence* pB) {
	try{
		
		vector<data_struct> data;
cout << q->getName() << '\t' << q->getAligned().length() << endl;		
		//vertical filter
		vector<Sequence*> temp;
		temp.push_back(q); temp.push_back(pA); temp.push_back(pB);
		verticalFilter(temp);
		
		//get these to avoid numerous function calls
		string query = q->getAligned();
		string parentA = pA->getAligned();
		string parentB = pB->getAligned();
		int length = query.length();
cout << q->getName() << '\t' << length << endl;
		
		//check window size
		if (length < (2*windowSize+windowStep)) { 
			mothurOut("Your window size is too large for " + q->getName() + ". I will make the window size " + toString(length/4) + " which is 1/4 the filtered length."); mothurOutEndLine();	
			windowSize = length / 4;
		}
		
		for (int i = windowSize-1; i <= (length - windowSize); i += windowStep) {
		
			int breakpoint = i;
			int leftLength = breakpoint + 1;
			int rightLength = length - leftLength;
			
			float QLA = computePercentID(query, parentA, 0, breakpoint);
			float QRB = computePercentID(query, parentB, breakpoint+1, length - 1);
		
			float QLB = computePercentID(query, parentB, 0, breakpoint);
			float QRA = computePercentID(query, parentA, breakpoint+1, length - 1);
		
			float LAB = computePercentID(parentA, parentB, 0, breakpoint);
			float RAB = computePercentID(parentA, parentB, breakpoint+1, length - 1);
		
			float AB = ((LAB*leftLength) + (RAB*rightLength)) / (float) length;
			float QA = ((QLA*leftLength) + (QRA*rightLength)) / (float) length;
			float QB = ((QLB*leftLength) + (QRB*rightLength)) / (float) length;
		
			float QLA_QRB = ((QLA*leftLength) + (QRB*rightLength)) / (float) length;
			float QLB_QRA = ((QLB*leftLength) + (QRA*rightLength)) / (float) length;
		
			//in original and not used
			//float avgQA_QB = ((QA*leftLength) + (QB*rightLength)) / (float) length;
		
			float divR_QLA_QRB = min((QLA_QRB/QA), (QLA_QRB/QB));
		
			float divR_QLB_QRA = min((QLB_QRA/QA), (QLB_QRA/QB));

			//is one of them above the 
			if (divR_QLA_QRB >= divRThreshold || divR_QLB_QRA >= divRThreshold) {
				
				if (((QLA_QRB > QA) && (QLA_QRB > QB) && (QLA >= parentFragmentThreshold) && (QRB >= parentFragmentThreshold))  ||
					((QLB_QRA > QA) && (QLB_QRA > QB) && (QLB >=parentFragmentThreshold) && (QRA >= parentFragmentThreshold)))  {
					
					data_struct member;
					
					member.divr_qla_qrb = divR_QLA_QRB;
					member.divr_qlb_qra = divR_QLB_QRA;
					member.qla_qrb = QLA_QRB;
					member.qlb_qra = QLB_QRA;
					member.qla = QLA;
					member.qrb = QRB;
					member.ab = AB; 
					member.qa = QA;
					member.qb = QB; 
					member.lab = LAB; 
					member.rab = RAB; 
					member.qra = QRA; 
					member.qlb = QLB; 
					member.winLStart = 0;
					member.winLEnd = breakpoint; 
					member.winRStart = breakpoint+1; 
					member.winREnd = length-1; 
					member.querySeq = *(q); 
					member.parentA = *(pA);
					member.parentB = *(pB);
					member.bsa = 0;
					member.bsb = 0;
					member.bsMax = 0;
					member.chimeraMax = 0;
					
					data.push_back(member);
					
				}//if
			}//if
		}//for
		
		return data;
		
	}
	catch(exception& e) {
		errorOut(e, "Slayer", "runBellerophon");
		exit(1);
	}
}
/***********************************************************************/
vector<snps> Slayer::getSNPS(string parentA, string query, string parentB, int left, int right) {
	try {
	
		vector<snps> data;

		for (int i = left; i <= right; i++) {
			
			char A = parentA[i];
			char Q = query[i];
			char B = parentB[i];
			
			if ((A != Q) || (B != Q)) {
				snps member;
				member.queryChar = Q;
				member.parentAChar = A;
				member.parentBChar = B;
				
				data.push_back(member);
			}
		}
		
		return data;
		
	}
	catch(exception& e) {
		errorOut(e, "Slayer", "getSNPS");
		exit(1);
	}
}
/***********************************************************************/
void Slayer::bootstrapSNPS(vector<snps> left, vector<snps> right, float& BSA, float& BSB) {
	try {

		srand((unsigned)time( NULL ));

		int count_A = 0; // sceneario QLA,QRB supported
		int count_B = 0; // sceneario QLB,QRA supported
	
		int numLeft = max(1, int(left.size()/10 +0.5));
		int numRight = max(1, int(right.size()/10 + 0.5));
		
		for (int i = 0; i < 100; i++) {
			//random sampling with replacement.
		
			vector<snps> selectedLeft;
			for (int j = 0; j < numLeft; j++) {
				int index = int(rand() % left.size());
				selectedLeft.push_back(left[index]);
			}
		
			vector<snps> selectedRight;
			for (int j = 0; j < numRight; j++) {
				int index = int(rand() % right.size());
				selectedRight.push_back(right[index]);
			}
		
			/* A  ------------------------------------------
			#       QLA                     QRA
			# Q  ------------------------------------------
			#                      |
			#                      |
			# Q  ------------------------------------------
			#       QLB                     QRB
			# B  ------------------------------------------ */
		
		
			float QLA = snpQA(selectedLeft);
			float QRA = snpQA(selectedRight);
		
			float QLB = snpQB(selectedLeft);
			float QRB = snpQB(selectedRight);
			
			//in original - not used - not sure why?
			//float ALB = snpAB(selectedLeft);
			//float ARB = snpAB(selectedRight);
		
			if ((QLA > QLB) && (QRB > QRA)) {
				count_A++;
			}
		
			if ((QLB > QLA) && (QRA > QRB)) {
				count_B++;
			}
		
		}

		BSA = (float) count_A;
		BSB = (float) count_B;
	
	}
	catch(exception& e) {
		errorOut(e, "Slayer", "bootstrapSNPS");
		exit(1);
	}
}
/***********************************************************************/
float Slayer::snpQA(vector<snps> data) {
	try {
	
		int numIdentical = 0;
	
		for (int i = 0; i < data.size(); i++) {
			if (data[i].parentAChar == data[i].queryChar) {
				numIdentical++;
			}
		}

		float percentID = (numIdentical / data.size()) * 100;
		
		return percentID;
	}
	catch(exception& e) {
		errorOut(e, "Slayer", "snpQA");
		exit(1);
	}
}
/***********************************************************************/
float Slayer::snpQB(vector<snps> data) {
	try {
	
		int numIdentical = 0;
	
		for (int i = 0; i < data.size(); i++) {
			if (data[i].parentBChar == data[i].queryChar) {
				numIdentical++;
			}
		}

		float percentID = (numIdentical / data.size()) * 100;
		
		return percentID;

	}
	catch(exception& e) {
		errorOut(e, "Slayer", "snpQB");
		exit(1);
	}
}
/***********************************************************************/
float Slayer::snpAB(vector<snps> data) {
	try {
		int numIdentical = 0;
	
		for (int i = 0; i < data.size(); i++) {
			if (data[i].parentAChar == data[i].parentBChar) {
				numIdentical++;
			}
		}

		float percentID = (numIdentical / data.size()) * 100;
		
		return percentID;

	}
	catch(exception& e) {
		errorOut(e, "Slayer", "snpAB");
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
//remove columns that contain any gaps
void Slayer::verticalFilter(vector<Sequence*> seqs) {
	try {
		vector<int> gaps;	gaps.resize(seqs[0]->getAligned().length(), 0);
		
		string filterString = (string(seqs[0]->getAligned().length(), '1'));
		
		//for each sequence
		for (int i = 0; i < seqs.size(); i++) {
		
			string seqAligned = seqs[i]->getAligned();
			
			for (int j = 0; j < seqAligned.length(); j++) {
				//if this spot is a gap
				if ((seqAligned[j] == '-') || (seqAligned[j] == '.') || (toupper(seqAligned[j]) == 'N'))	{	gaps[j]++;	}
			}
		}
		
		//zero out spot where all sequences have blanks
		int numColRemoved = 0;

		for(int i = 0; i < seqs[0]->getAligned().length(); i++){
			if(gaps[i] != 0)	{	filterString[i] = '0'; 	numColRemoved++;  }
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
