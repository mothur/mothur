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
Slayer::Slayer(int win, int increment, int parentThreshold, float div, int i, int snp, int mi) :
		minBS(mi), windowSize(win), windowStep(increment), parentFragmentThreshold(parentThreshold), divRThreshold(div), iters(i), percentSNPSample(snp){ m = MothurOut::getInstance(); }
/***********************************************************************/
string Slayer::getResults(Sequence query, vector<Sequence> refSeqs) {
	try {
		vector<data_struct> all; all.clear();
		myQuery = query;

		for (int i = 0; i < refSeqs.size(); i++) {
		
			for (int j = i+1; j < refSeqs.size(); j++) {
			
				if (m->getControl_pressed()) { return "no";  }
	
				//make copies of query and each parent because runBellerophon removes gaps and messes them up
				Sequence q(query.getName(), query.getAligned());
				Sequence leftParent(refSeqs[i].getName(), refSeqs[i].getAligned());
				Sequence rightParent(refSeqs[j].getName(), refSeqs[j].getAligned());
				
	
				map<int, int> spots;  //map from spot in original sequence to spot in filtered sequence for query and both parents
				vector<data_struct> divs = runBellerophon(q, leftParent, rightParent, spots);
	
				if (m->getControl_pressed()) { return "no"; }

				vector<data_struct> selectedDivs;
				for (int k = 0; k < divs.size(); k++) {
					
					vector<snps> snpsLeft = getSNPS(divs[k].parentA.getAligned(), divs[k].querySeq.getAligned(), divs[k].parentB.getAligned(), divs[k].winLStart, divs[k].winLEnd);
					vector<snps> snpsRight = getSNPS(divs[k].parentA.getAligned(), divs[k].querySeq.getAligned(), divs[k].parentB.getAligned(), divs[k].winRStart, divs[k].winREnd);
	
					if (m->getControl_pressed()) {  return "no"; }
					
					int numSNPSLeft = snpsLeft.size();
					int numSNPSRight = snpsRight.size();
					
					//require at least 4 SNPs on each side of the break
					if ((numSNPSLeft >= 4) && (numSNPSRight >= 4)) {
							
						float BS_A, BS_B;
						bootstrapSNPS(snpsLeft, snpsRight, BS_A, BS_B, iters);
						
						if (m->getControl_pressed()) { return "no"; }

						divs[k].bsa = BS_A;
						divs[k].bsb = BS_B;
						divs[k].bsMax = max(BS_A, BS_B);
						divs[k].chimeraMax = max(divs[k].qla_qrb, divs[k].qlb_qra);
						
						
						//so results reflect orignal alignment
						divs[k].winLStart = spots[divs[k].winLStart];
						divs[k].winLEnd = spots[divs[k].winLEnd];  
						divs[k].winRStart = spots[divs[k].winRStart]; 
						divs[k].winREnd = spots[divs[k].winREnd]; 
						
						selectedDivs.push_back(divs[k]);
					}
				}

				//save selected
				for (int mi = 0; mi < selectedDivs.size(); mi++) {  all.push_back(selectedDivs[mi]);	}
			}
		}
		

		// compute bootstrap support
		if (all.size() > 0) {
			//sort them
			sort(all.begin(), all.end(), compareDataStruct);
			reverse(all.begin(), all.end());
						
			outputResults = all;
			return "yes"; 
		}else {
			outputResults = all;
			return "no";
		}
	}
	catch(exception& e) {
		m->errorOut(e, "Slayer", "getResults");
		exit(1);
	}
}
/***********************************************************************/
vector<data_struct> Slayer::runBellerophon(Sequence q, Sequence pA, Sequence pB, map<int, int>& spots) {
	try{
		
		vector<data_struct> data;
				
		//maps spot in new alignment to spot in alignment before filter
		spots = verticalFilter(q, pA, pB);  //fills baseSpots
		
		//get these to avoid numerous function calls
		string query = q.getAligned();
		string parentA = pA.getAligned();
		string parentB = pB.getAligned();
		int length = query.length();
	
		//check window size
		if (length < (2*windowSize+windowStep)) { 
//			m->mothurOut("Your window size is too large for " + q->getName() + ". I will make the window size " + toString(length/4) + " which is 1/4 the filtered length.\n"); 	
			windowSize = length / 4;
		}
		
		for (int i = windowSize-1; i <= (length - windowSize); i += windowStep) {
		
			if (m->getControl_pressed()) { return data; }
		
			int breakpoint = i;
			int leftLength = breakpoint + 1;
			int rightLength = length - leftLength;
				
			float QLA = computePercentID(query, parentA, 0, breakpoint);
			float QRB = computePercentID(query, parentB, breakpoint+1, length-1);
		
			float QLB = computePercentID(query, parentB, 0, breakpoint);
			float QRA = computePercentID(query, parentA, breakpoint+1, length-1);
		
			float LAB = computePercentID(parentA, parentB, 0, breakpoint);
			float RAB = computePercentID(parentA, parentB, breakpoint+1, length-1);	
			
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
					member.querySeq = q;
					member.parentA = pA;
					member.parentB = pB;
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
		m->errorOut(e, "Slayer", "runBellerophon");
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

				//ensure not neighboring a gap. change to 12/09 release of chimeraSlayer - not sure what this adds, but it eliminates alot of SNPS

				
				if (
					//did query loose a base here during filter??
					( i == 0 || abs (baseSpots[0][i] - baseSpots[0][i-1]) == 1) &&
					( i == query.length()-1 || abs (baseSpots[0][i] - baseSpots[0][i+1]) == 1)
					&&
					//did parentA loose a base here during filter??
					( i == 0 || abs (baseSpots[1][i] - baseSpots[1][i-1]) == 1) &&
					( i == parentA.length()-1 || abs (baseSpots[1][i] - baseSpots[1][i+1]) == 1) 
					&&
					//did parentB loose a base here during filter??
					( i == 0 || abs (baseSpots[2][i] - baseSpots[2][i-1]) == 1) &&
					( i == parentB.length()-1 || abs (baseSpots[2][i] - baseSpots[2][i+1]) == 1)
					) 
				{ 
					snps member;
					member.queryChar = Q;
					member.parentAChar = A;
					member.parentBChar = B;
					data.push_back(member);
				}
			}
		}
		
		return data;
		
	}
	catch(exception& e) {
		m->errorOut(e, "Slayer", "getSNPS");
		exit(1);
	}
}
/***********************************************************************/
int Slayer::bootstrapSNPS(vector<snps> left, vector<snps> right, float& BSA, float& BSB, int numIters) {
	try {

		m->setRandomSeed((unsigned)time( nullptr ));

		int count_A = 0; // sceneario QLA,QRB supported
		int count_B = 0; // sceneario QLB,QRA supported
	
		int numLeft = max(1, int(left.size() * percentSNPSample/(float)100 + 0.5));
		int numRight = max(1, int(right.size() * percentSNPSample/(float)100 + 0.5));
        Utils util;
		for (int i = 0; i < numIters; i++) {
			//random sampling with replacement.
		
			if (m->getControl_pressed()) { return 0;  }
			
			vector<snps> selectedLeft;

			for (int j = 0; j < numLeft; j++) {
				int index = util.getRandomIndex((int)left.size()-1);
				selectedLeft.push_back(left[index]);
			}

			vector<snps> selectedRight;
			for (int j = 0; j < numRight; j++) {
				int index = util.getRandomIndex((int)right.size()-1);
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

		BSA = (float) count_A / (float) numIters * 100;
		BSB = (float) count_B / (float) numIters * 100;

		return 0;
	
	}
	catch(exception& e) {
		m->errorOut(e, "Slayer", "bootstrapSNPS");
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

		float percentID = (numIdentical / (float) data.size()) * 100;
		
		return percentID;
	}
	catch(exception& e) {
		m->errorOut(e, "Slayer", "snpQA");
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

		float percentID = (numIdentical / (float) data.size()) * 100;
		
		return percentID;

	}
	catch(exception& e) {
		m->errorOut(e, "Slayer", "snpQB");
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

		float percentID = (numIdentical / (float) data.size()) * 100;
		
		return percentID;

	}
	catch(exception& e) {
		m->errorOut(e, "Slayer", "snpAB");
		exit(1);
	}
}
/***********************************************************************/
float Slayer::computePercentID(string queryAlign, string chimera, int left, int right) {
	try {
				
		int numIdentical = 0;
		int countA = 0;
		int countB = 0;
		for (int i = left; i <= right; i++) {
			if (((queryAlign[i] != 'G') && (queryAlign[i] != 'T') && (queryAlign[i] != 'A') && (queryAlign[i] != 'C')&& (queryAlign[i] != '.') && (queryAlign[i] != '-')) ||
				((chimera[i] != 'G') && (chimera[i] != 'T') && (chimera[i] != 'A') && (chimera[i] != 'C')&& (chimera[i] != '.') && (chimera[i] != '-'))) {}
			else {
				
				bool charA = false; bool charB = false;
				if ((queryAlign[i] == 'G') || (queryAlign[i] == 'T') || (queryAlign[i] == 'A') || (queryAlign[i] == 'C')) { charA = true; }
				if ((chimera[i] == 'G') || (chimera[i] == 'T') || (chimera[i] == 'A') || (chimera[i] == 'C')) { charB = true; }
				
				if (charA || charB) {
					
					if (charA) { countA++; }
					if (charB) { countB++; }
					
					if (queryAlign[i] == chimera[i]) {
						numIdentical++;
					}
				}
			}
			
		}
		
		float numBases = (countA + countB) /(float) 2;
		
		if (numBases == 0) { return 0; }
		
		float percentIdentical = (numIdentical/(float)numBases) * 100;

		return percentIdentical;
		
	}
	catch(exception& e) {
		m->errorOut(e, "Slayer", "computePercentID");
		exit(1);
	}
}
/***********************************************************************/
//remove columns that contain any gaps
map<int, int> Slayer::verticalFilter(Sequence& q, Sequence& pA, Sequence& pB) {
	try {
		//find baseSpots
		baseSpots.clear(); 
		baseSpots.resize(3);  //query, parentA, parentB
	
		vector<int> gaps;	gaps.resize(q.getAligned().length(), 0);
		
		string filterString = (string(q.getAligned().length(), '1'));
		
		string seqAligned = q.getAligned();
		for (int j = 0; j < seqAligned.length(); j++) {
			//if this spot is a gap
			if ((seqAligned[j] == '-') || (seqAligned[j] == '.') || (toupper(seqAligned[j]) == 'N'))	{   gaps[j]++;	}
		}
		
		seqAligned = pA.getAligned();
		for (int j = 0; j < seqAligned.length(); j++) {
			//if this spot is a gap
			if ((seqAligned[j] == '-') || (seqAligned[j] == '.') || (toupper(seqAligned[j]) == 'N'))	{   gaps[j]++;	}
		}
		
		seqAligned = pB.getAligned();
		for (int j = 0; j < seqAligned.length(); j++) {
			//if this spot is a gap
			if ((seqAligned[j] == '-') || (seqAligned[j] == '.') || (toupper(seqAligned[j]) == 'N'))	{   gaps[j]++;	}
		}
		
		
		//zero out spot where any sequences have blanks
		int numColRemoved = 0;
		int count = 0;
		map<int, int> maskMap; maskMap.clear();

		for(int i = 0; i < q.getAligned().length(); i++){
			if(gaps[i] != 0)	{	filterString[i] = '0'; 	numColRemoved++;  }
			else {
				maskMap[count] = i;
				count++;
			}
		}

		seqAligned = q.getAligned();
		string newAligned = "";
			
		int baseCount = 0;
		count = 0;
		for (int j = 0; j < seqAligned.length(); j++) {
			//are you a base
			if ((seqAligned[j] != '-') && (seqAligned[j] != '.') && (toupper(seqAligned[j]) != 'N'))	{ baseCount++; }
			
			//if this spot is not a gap
			if (filterString[j] == '1') { 
				newAligned += seqAligned[j]; 
				baseSpots[0][count] = baseCount;
				count++;
			}
		}
			
		q.setAligned(newAligned);
		
		seqAligned = pA.getAligned();
		newAligned = "";
		
		baseCount = 0;
		count = 0;
		for (int j = 0; j < seqAligned.length(); j++) {
			//are you a base
			if ((seqAligned[j] != '-') && (seqAligned[j] != '.') && (toupper(seqAligned[j]) != 'N'))	{ baseCount++; }
			
			//if this spot is not a gap
			if (filterString[j] == '1') { 
				newAligned += seqAligned[j]; 
				baseSpots[1][count] = baseCount;
				count++;
			}
		}
		
		pA.setAligned(newAligned);
		
		seqAligned = pB.getAligned();
		newAligned = "";
		
		baseCount = 0;
		count = 0;
		for (int j = 0; j < seqAligned.length(); j++) {
			//are you a base
			if ((seqAligned[j] != '-') && (seqAligned[j] != '.') && (toupper(seqAligned[j]) != 'N'))	{ baseCount++; }
			
			//if this spot is not a gap
			if (filterString[j] == '1') { 
				newAligned += seqAligned[j]; 
				baseSpots[2][count] = baseCount;
				count++;
			}
		}
		
		pB.setAligned(newAligned);
		
		
		return maskMap;
	}
	catch(exception& e) {
		m->errorOut(e, "Slayer", "verticalFilter");
		exit(1);
	}
}
/***********************************************************************/
