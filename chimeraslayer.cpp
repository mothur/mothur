/*
 *  chimeraslayer.cpp
 *  Mothur
 *
 *  Created by westcott on 9/25/09.
 *  Copyright 2009 Pschloss Lab. All rights reserved.
 *
 */

#include "chimeraslayer.h"
#include "chimerarealigner.h"

//***************************************************************************************************************
ChimeraSlayer::ChimeraSlayer(string mode) : searchMethod(mode) {  	decalc = new DeCalculator();	  }
//***************************************************************************************************************
ChimeraSlayer::~ChimeraSlayer() { 	delete decalc;	 }
//***************************************************************************************************************
void ChimeraSlayer::printHeader(ostream& out) {
	mothurOutEndLine();
	mothurOut("Only reporting sequence supported by 90% of bootstrapped results.");
	mothurOutEndLine();
}
//***************************************************************************************************************
void ChimeraSlayer::print(ostream& out) {
	try {
		if (chimeraFlags == "yes") {
			string chimeraFlag = "no";
			if(  (chimeraResults[0].bsa >= minBS && chimeraResults[0].divr_qla_qrb >= divR)
			   ||
			   (chimeraResults[0].bsb >= minBS && chimeraResults[0].divr_qlb_qra >= divR) ) { chimeraFlag = "yes"; }
			
			
			if (chimeraFlag == "yes") {	
				if ((chimeraResults[0].bsa >= minBS) || (chimeraResults[0].bsb >= minBS)) {
					mothurOut(querySeq->getName() + "\tyes"); mothurOutEndLine();
				}
			}
			out << querySeq->getName() << "\tyes" << endl;
			printBlock(chimeraResults[0], out);
			out << endl;
		}else {  out << querySeq->getName() << "\tno" << endl;  }
		
	}
	catch(exception& e) {
		errorOut(e, "ChimeraSlayer", "print");
		exit(1);
	}
}
//***************************************************************************************************************
int ChimeraSlayer::getChimeras(Sequence* query) {
	try {
		chimeraFlags = "no";
		querySeq = query;
		
		for (int i = 0; i < query->getAligned().length(); i++) {
			spotMap[i] = i;
		}
		
		//referenceSeqs, numWanted, matchScore, misMatchPenalty, divR, minSimilarity
		maligner = new Maligner(templateSeqs, numWanted, match, misMatch, divR, minSim, minCov, searchMethod);
		slayer = new Slayer(window, increment, minSim, divR, iters, minSNP);
		
		string chimeraFlag = maligner->getResults(query, decalc);
		vector<results> Results = maligner->getOutput();

		//realign query to parents to improve slayers detection rate
		ChimeraReAligner realigner(templateSeqs, match, misMatch);
		realigner.reAlign(query, Results);
cout << query->getName() << '\n' << query->getAligned() << endl;
			//if (chimeraFlag == "yes") {
			
		//get sequence that were given from maligner results
		vector<SeqDist> seqs;
		map<string, float> removeDups;
		map<string, float>::iterator itDup;
		for (int j = 0; j < Results.size(); j++) {
			float dist = (Results[j].regionEnd - Results[j].regionStart + 1) * Results[j].queryToParentLocal;
			//only add if you are not a duplicate
			itDup = removeDups.find(Results[j].parent);
			if (itDup == removeDups.end()) { //this is not duplicate
				removeDups[Results[j].parent] = dist;
			}else if (dist > itDup->second) { //is this a stronger number for this parent
				removeDups[Results[j].parent] = dist;
			}
		}
		
		for (itDup = removeDups.begin(); itDup != removeDups.end(); itDup++) {
			Sequence* seq = getSequence(itDup->first); //makes copy so you can filter and mask and not effect template
			
			SeqDist member;
			member.seq = seq;
			member.dist = itDup->second;
			
			seqs.push_back(member);
		}
		
		//limit number of parents to explore - default 3
		if (Results.size() > parents) {
			//sort by distance
			sort(seqs.begin(), seqs.end(), compareSeqDist);
			//prioritize larger more similiar sequence fragments
			reverse(seqs.begin(), seqs.end());
			
			for (int k = seqs.size()-1; k > (parents-1); k--)  {  
				delete seqs[k].seq;
				seqs.pop_back();	
			}
		}
		
		//put seqs into vector to send to slayer
		vector<Sequence*> seqsForSlayer;
		for (int k = 0; k < seqs.size(); k++) {  seqsForSlayer.push_back(seqs[k].seq);	}
		//cout << i+1 << "num parents = " << seqsForSlayer.size() << '\t' << chimeraFlag << endl;
//ofstream out;
//string name = querySeqs[i]->getName();
//cout << name << endl;
//name = name.substr(name.find_first_of("|")+1);
//cout << name << endl;
//name = name.substr(name.find_first_of("|")+1);
//cout << name << endl;
//name = name.substr(0, name.find_last_of("|"));
//cout << name << endl;
//string filename = toString(i+1) + ".seqsforslayer";
//openOutputFile(filename, out);	
//cout << querySeqs[i]->getName() << endl;
//for (int u = 0; u < seqsForSlayer.size(); u++) { cout << seqsForSlayer[u]->getName() << '\t'; seqsForSlayer[u]->printSequence(out);	}
//cout << endl;
//out.close();
//filename = toString(i+1) + ".fasta";
//openOutputFile(filename, out);	
//querySeqs[i]->printSequence(out);
//out.close();


			
		//mask then send to slayer...
		if (seqMask != "") {
			decalc->setMask(seqMask);
			
			//mask querys
			decalc->runMask(query);
			
			//mask parents
			for (int k = 0; k < seqsForSlayer.size(); k++) {
				decalc->runMask(seqsForSlayer[k]);
			}
			
			spotMap = decalc->getMaskMap();
		}
		
		//send to slayer
		chimeraFlags = slayer->getResults(query, seqsForSlayer);
		chimeraResults = slayer->getOutput();
	
		//free memory
		for (int k = 0; k < seqs.size(); k++) {  delete seqs[k].seq;   }
		//}
			
		return 0;
	}
	catch(exception& e) {
		errorOut(e, "ChimeraSlayer", "getChimeras");
		exit(1);
	}
}
//***************************************************************************************************************
void ChimeraSlayer::printBlock(data_struct data, ostream& out){
	try {
		
		out << "parentA: " << data.parentA.getName() << "  parentB: " << data.parentB.getName()  << endl;
		out << "Left Window: " << spotMap[data.winLStart] << " " << spotMap[data.winLEnd] << endl;
		out << "Right Window: " << spotMap[data.winRStart] << " " << spotMap[data.winREnd] << endl;
		
		out << "Divergence of Query to Leftside ParentA and Rightside ParentB: " << data.divr_qla_qrb << '\t' << "Bootstrap: " << data.bsa << endl;
		out << "Divergence of Query to Rightside ParentA and Leftside ParentB: " << data.divr_qlb_qra << '\t' << "Bootstrap: " << data.bsb << endl;
		
		out << "Similarity of parents: " << data.ab << endl;
		out << "Similarity of query to parentA: " << data.qa << endl;
		out << "Similarity of query to parentB: " << data.qb << endl;
		
		out << "Percent_ID QLA_QRB: " << data.qla_qrb << endl;
		out << "Percent_ID QLB_QRA: " << data.qlb_qra << endl;
		//out << "Per_id(QL,A): " << data.qla << endl;
		//out << "Per_id(QL,B): " << data.qlb << endl;
		//out << "Per_id(QR,A): " << data.qra << endl;
		//out << "Per_id(QR,B): " << data.qrb << endl;

		
		out << "DeltaL: " << (data.qla - data.qlb) << endl;
		out << "DeltaR: " << (data.qra - data.qrb) << endl;

	}
	catch(exception& e) {
		errorOut(e, "ChimeraSlayer", "printBlock");
		exit(1);
	}
}
//***************************************************************************************************************

