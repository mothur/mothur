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
string Slayer::getResults(Sequence* query, vector<Sequence*> refSeqs) {
	try {
		vector<data_struct> all; all.clear();
		myQuery = *query;
		
		for (int i = 0; i < refSeqs.size(); i++) {
		
			for (int j = i+1; j < refSeqs.size(); j++) {
			
				if (m->control_pressed) { return "no";  }
	
				//make copies of query and each parent because runBellerophon removes gaps and messes them up
				Sequence* q = new Sequence(query->getName(), query->getAligned());
				Sequence* leftParent = new Sequence(refSeqs[i]->getName(), refSeqs[i]->getAligned());
				Sequence* rightParent = new Sequence(refSeqs[j]->getName(), refSeqs[j]->getAligned());
		//cout << "parents: (" << refSeqs[i]->getName() << ", " << refSeqs[j]->getName() << ")\n";
				map<int, int> spots;  //map from spot in original sequence to spot in filtered sequence for query and both parents
				vector<data_struct> divs = runBellerophon(q, leftParent, rightParent, spots);
				//cout << divs.size() << endl;
				if (m->control_pressed) { 
					delete q;
					delete leftParent;
					delete rightParent;
					return "no"; 
				}
					
//				cout << divs.size() << endl;
				vector<data_struct> selectedDivs;
				for (int k = 0; k < divs.size(); k++) {
					
					vector<snps> snpsLeft = getSNPS(divs[k].parentA.getAligned(), divs[k].querySeq.getAligned(), divs[k].parentB.getAligned(), divs[k].winLStart, divs[k].winLEnd);
					vector<snps> snpsRight = getSNPS(divs[k].parentA.getAligned(), divs[k].querySeq.getAligned(), divs[k].parentB.getAligned(), divs[k].winRStart, divs[k].winREnd);
					//cout << refSeqs[i]->getName() << '\t' << refSeqs[j]->getName() << '\t' << k << divs[k].parentA.getAligned() << endl << divs[k].parentB.getAligned() << endl;	
					if (m->control_pressed) { delete q; delete leftParent; delete rightParent; return "no"; }
					
					int numSNPSLeft = snpsLeft.size();
					int numSNPSRight = snpsRight.size();
					
					//require at least 4 SNPs on each side of the break
					if ((numSNPSLeft >= 4) && (numSNPSRight >= 4)) {
							
						float BS_A, BS_B;
						bootstrapSNPS(snpsLeft, snpsRight, BS_A, BS_B, iters);
						
						if (m->control_pressed) { delete q; delete leftParent; delete rightParent; return "no"; }

						divs[k].bsa = BS_A;
						divs[k].bsb = BS_B;
						divs[k].bsMax = max(BS_A, BS_B);
						divs[k].chimeraMax = max(divs[k].qla_qrb, divs[k].qlb_qra);
						
						
						//are we within 10 points of the bootstrap cutoff?
						if ((divs[k].bsMax >= (minBS-10)) && (iters < 1000)) {
							bootstrapSNPS(snpsLeft, snpsRight, BS_A, BS_B, 1000);
								
							if (m->control_pressed) { delete q; delete leftParent; delete rightParent; return "no"; }
								
							divs[k].bsa = BS_A;
							divs[k].bsb = BS_B;
							divs[k].bsMax = max(BS_A, BS_B);
							divs[k].chimeraMax = max(divs[k].qla_qrb, divs[k].qlb_qra);
						}
						
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
vector<data_struct> Slayer::runBellerophon(Sequence* q, Sequence* pA, Sequence* pB, map<int, int>& spots) {
	try{
		
		vector<data_struct> data;
		
		//vertical filter
		vector<Sequence*> temp;
		temp.push_back(q); temp.push_back(pA); temp.push_back(pB);
		
		//maps spot in new alignment to spot in alignment before filter
		spots = verticalFilter(temp);  //fills baseSpots
		
		//get these to avoid numerous function calls
		string query = q->getAligned();
		string parentA = pA->getAligned();
		string parentB = pB->getAligned();
		int length = query.length();
//cout << q->getName() << endl << q->getAligned() << endl << endl;	
//cout << pA->getName() << endl << pA->getUnaligned() << endl << endl;		
//cout << pB->getName() << endl << pB->getUnaligned() << endl << endl;	
//cout << " length = " << length << endl;
	
		//check window size
		if (length < (2*windowSize+windowStep)) { 
//			m->mothurOut("Your window size is too large for " + q->getName() + ". I will make the window size " + toString(length/4) + " which is 1/4 the filtered length."); m->mothurOutEndLine();	
			windowSize = length / 4;
		}
		
		for (int i = windowSize-1; i <= (length - windowSize); i += windowStep) {
		
			if (m->control_pressed) { return data; }
		
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
			//cout << q->getName() << '\t';
			//cout << pA->getName() << '\t';
			//cout << pB->getName() << '\t';
		   // cout << "bp: " << breakpoint << " CHIM_TYPE_A\t" << divR_QLA_QRB << "\tQLA: " << QLA << "\tQRB: " << QRB << "\tQLA_QRB: " << QLA_QRB;
			//cout << "\tCHIM_TYPE_B\t" << divR_QLB_QRA << "\tQLB: " << QLB << "\tQRA: " << QRA << "\tQLB_QRA: " << QLB_QRA << endl;
//cout << leftLength << '\t' << rightLength << '\t' << QLA << '\t' << QRB << '\t' << QLB << '\t' << QRA  << '\t' << LAB << '\t' << RAB << '\t' << AB << '\t' << QA << '\t' << QB << '\t' << QLA_QRB << '\t' <<  QLB_QRA <<    endl;    		

//cout << divRThreshold << endl;
//cout << breakpoint << '\t' << divR_QLA_QRB << '\t' << divR_QLB_QRA << endl;
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
		m->errorOut(e, "Slayer", "runBellerophon");
		exit(1);
	}
}
/***********************************************************************/
vector<snps> Slayer::getSNPS(string parentA, string query, string parentB, int left, int right) {
	try {
	
		vector<snps> data;
//cout << left << '\t' << right << endl;
		for (int i = left; i <= right; i++) {
			
			char A = parentA[i];
			char Q = query[i];
			char B = parentB[i];
			
			if ((A != Q) || (B != Q)) {
//cout << "not equal " << Q << '\t' << A << '\t' << B << endl;
			
				//ensure not neighboring a gap. change to 12/09 release of chimeraSlayer - not sure what this adds, but it eliminates alot of SNPS
				if (
					//did query loose a base here during filter??
					( i == 0 || abs (baseSpots[0][i] - baseSpots[0][i-1]) == 1) &&
					( i == query.length() || abs (baseSpots[0][i] - baseSpots[0][i+1]) == 1)
					&&
					//did parentA loose a base here during filter??
					( i == 0 || abs (baseSpots[1][i] - baseSpots[1][i-1]) == 1) &&
					( i == parentA.length() || abs (baseSpots[1][i] - baseSpots[1][i+1]) == 1) 
					&&
					//did parentB loose a base here during filter??
					( i == 0 || abs (baseSpots[2][i] - baseSpots[2][i-1]) == 1) &&
					( i == parentB.length() || abs (baseSpots[2][i] - baseSpots[2][i+1]) == 1)
					) 
				{ 
				
					snps member;
					member.queryChar = Q;
					member.parentAChar = A;
					member.parentBChar = B;
//cout << "not neighboring a gap " << Q << '\t' << A << '\t' << B << '\t' << baseSpots[0][i] << '\t' << baseSpots[0][i+1] << '\t' << baseSpots[0][i-1] << '\t' << baseSpots[1][i] << '\t' << baseSpots[1][i+1] << '\t' << baseSpots[1][i-1] << '\t' << baseSpots[2][i] << '\t' << baseSpots[2][i+1] << '\t' << baseSpots[2][i-1] << endl;				
					data.push_back(member);
				}
			}
//			cout << i << '\t' << data.size() << endl;
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

		srand((unsigned)time( NULL ));

		int count_A = 0; // sceneario QLA,QRB supported
		int count_B = 0; // sceneario QLB,QRA supported
	
		int numLeft = max(1, int(left.size() * percentSNPSample/(float)100 + 0.5));
		int numRight = max(1, int(right.size() * percentSNPSample/(float)100 + 0.5));
		//cout << numLeft << '\t' << numRight << endl;
		for (int i = 0; i < numIters; i++) {
			//random sampling with replacement.
		
			if (m->control_pressed) { return 0;  }
			
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
			
//cout << "selected left snp: \n";
//for (int j = 0; j < selectedLeft.size(); j++) {  cout << selectedLeft[j].parentAChar;  } 
//cout << endl;
//for (int j = 0; j < selectedLeft.size(); j++) {  cout << selectedLeft[j].queryChar;  }
//cout << endl;
//for (int j = 0; j < selectedLeft.size(); j++) {  cout << selectedLeft[j].parentBChar;  }
//cout << endl;
//cout << "selected right snp: \n";
//for (int j = 0; j < selectedRight.size(); j++) {  cout << selectedRight[j].parentAChar;  } 
//cout << endl;
//for (int i = 0; i < selectedRight.size(); i++) {  cout << selectedRight[i].queryChar;  }
//cout << endl;
//for (int i = 0; i < selectedRight.size(); i++) {  cout << selectedRight[i].parentBChar;  }
//cout << endl;		
		}


		//cout << count_A << '\t' << count_B << endl;

		BSA = (float) count_A / (float) numIters * 100;
		BSB = (float) count_B / (float) numIters * 100;
//cout << "bsa = " << BSA << " bsb = " << BSB << endl;

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
		m->errorOut(e, "Slayer", "computePercentID");
		exit(1);
	}
}
/***********************************************************************/
//remove columns that contain any gaps
map<int, int> Slayer::verticalFilter(vector<Sequence*> seqs) {
	try {
		//find baseSpots
		baseSpots.clear(); 
		baseSpots.resize(3);  //query, parentA, parentB
	
		vector<int> gaps;	gaps.resize(seqs[0]->getAligned().length(), 0);
		
		string filterString = (string(seqs[0]->getAligned().length(), '1'));
		
		//for each sequence
		for (int i = 0; i < seqs.size(); i++) {
		
			string seqAligned = seqs[i]->getAligned();
			
			for (int j = 0; j < seqAligned.length(); j++) {
				//if this spot is a gap
				if ((seqAligned[j] == '-') || (seqAligned[j] == '.') || (toupper(seqAligned[j]) == 'N'))	{   gaps[j]++;	}
			}
		}
		
		//zero out spot where any sequences have blanks
		int numColRemoved = 0;
		int count = 0;
		map<int, int> maskMap; maskMap.clear();

		for(int i = 0; i < seqs[0]->getAligned().length(); i++){
			if(gaps[i] != 0)	{	filterString[i] = '0'; 	numColRemoved++;  }
			else {
				maskMap[count] = i;
				count++;
			}
		}

		//for each sequence
		for (int i = 0; i < seqs.size(); i++) {
		
			string seqAligned = seqs[i]->getAligned();
			string newAligned = "";
			
			int baseCount = 0;
			int count = 0;
			for (int j = 0; j < seqAligned.length(); j++) {
				//are you a base
				if ((seqAligned[j] != '-') && (seqAligned[j] != '.') && (toupper(seqAligned[j]) != 'N'))	{ baseCount++; }
			
				//if this spot is not a gap
				if (filterString[j] == '1') { 
					newAligned += seqAligned[j]; 
					baseSpots[i][count] = baseCount;
					count++;
				}
			}
			
			seqs[i]->setAligned(newAligned);
		}
		
		return maskMap;
	}
	catch(exception& e) {
		m->errorOut(e, "Slayer", "verticalFilter");
		exit(1);
	}
}
/***********************************************************************/
