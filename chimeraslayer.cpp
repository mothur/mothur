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
#include "kmerdb.hpp"

//***************************************************************************************************************
ChimeraSlayer::ChimeraSlayer(string mode, bool r, string f) : searchMethod(mode), realign(r), fastafile(f) {  	
	decalc = new DeCalculator();	
}
//***************************************************************************************************************
int ChimeraSlayer::doPrep() {
	try {
	
		string 	kmerDBNameLeft;
		string 	kmerDBNameRight;
		
		//generate the kmerdb to pass to maligner
		if (searchMethod == "kmer") { 
			//leftside
			string leftTemplateFileName = "left." + templateFileName;
			databaseLeft = new KmerDB(leftTemplateFileName, kmerSize);			
			kmerDBNameLeft = leftTemplateFileName.substr(0,leftTemplateFileName.find_last_of(".")+1) + char('0'+ kmerSize) + "mer";
			ifstream kmerFileTestLeft(kmerDBNameLeft.c_str());
			
			if(!kmerFileTestLeft){	
			
				for (int i = 0; i < templateSeqs.size(); i++) {
					
					if (m->control_pressed) { return 0; } 
					
					string leftFrag = templateSeqs[i]->getUnaligned();
					leftFrag = leftFrag.substr(0, int(leftFrag.length() * 0.33));
					
					Sequence leftTemp(templateSeqs[i]->getName(), leftFrag);
					databaseLeft->addSequence(leftTemp);	
				}
				databaseLeft->generateDB();
				
			}else {	
				databaseLeft->readKmerDB(kmerFileTestLeft);	
			}
			kmerFileTestLeft.close();
			
			databaseLeft->setNumSeqs(templateSeqs.size());
			
			//rightside
			string rightTemplateFileName = "right." + templateFileName;
			databaseRight = new KmerDB(rightTemplateFileName, kmerSize);			
			kmerDBNameRight = rightTemplateFileName.substr(0,rightTemplateFileName.find_last_of(".")+1) + char('0'+ kmerSize) + "mer";
			ifstream kmerFileTestRight(kmerDBNameRight.c_str());
			
			if(!kmerFileTestRight){	
			
				for (int i = 0; i < templateSeqs.size(); i++) {
					if (m->control_pressed) { return 0; } 
					
					string rightFrag = templateSeqs[i]->getUnaligned();
					rightFrag = rightFrag.substr(int(rightFrag.length() * 0.66));
					
					Sequence rightTemp(templateSeqs[i]->getName(), rightFrag);
					databaseRight->addSequence(rightTemp);	
				}
				databaseRight->generateDB();
				
			}else {	
				databaseRight->readKmerDB(kmerFileTestRight);	
			}
			kmerFileTestRight.close();
			
			databaseRight->setNumSeqs(templateSeqs.size());

		}
		
		int start = time(NULL);	
		//filter the sequences
		//read in all query seqs
		ifstream in; 
		openInputFile(fastafile, in);
		
		vector<Sequence*> tempQuerySeqs;
		while(!in.eof()){
			if (m->control_pressed) { for (int i = 0; i < tempQuerySeqs.size(); i++) { delete tempQuerySeqs[i];  } return 0; } 
		
			Sequence* s = new Sequence(in);
			gobble(in);
			
			if (s->getName() != "") { tempQuerySeqs.push_back(s); }
		}
		in.close();
		
		vector<Sequence*> temp = templateSeqs;
		for (int i = 0; i < tempQuerySeqs.size(); i++) {  temp.push_back(tempQuerySeqs[i]);  }
				
		createFilter(temp, 0.0); //just removed columns where all seqs have a gap
				
		for (int i = 0; i < tempQuerySeqs.size(); i++) { delete tempQuerySeqs[i];  }
		
		if (m->control_pressed) {  return 0; } 

		
		//run filter on template
		for (int i = 0; i < templateSeqs.size(); i++) {  if (m->control_pressed) {  return 0; }  runFilter(templateSeqs[i]);  }
		
		m->mothurOutEndLine(); m->mothurOut("It took " + toString(time(NULL) - start) + " secs to filter.");	m->mothurOutEndLine();
		
		return 0;

	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayer", "doprep");
		exit(1);
	}
}
//***************************************************************************************************************
ChimeraSlayer::~ChimeraSlayer() { 	delete decalc;  if (searchMethod == "kmer") {  delete databaseRight;  delete databaseLeft;  }	 }
//***************************************************************************************************************
void ChimeraSlayer::printHeader(ostream& out) {
	m->mothurOutEndLine();
	m->mothurOut("Only reporting sequence supported by " + toString(minBS) + "% of bootstrapped results.");
	m->mothurOutEndLine();
	
	out << "Name\tLeftParent\tRightParent\tDivQLAQRB\tPerIDQLAQRB\tBootStrapA\tDivQLBQRA\tPerIDQLBQRA\tBootStrapB\tFlag\tLeftWindow\tRightWindow\n";
}
//***************************************************************************************************************
int ChimeraSlayer::print(ostream& out, ostream& outAcc) {
	try {
		if (chimeraFlags == "yes") {
			string chimeraFlag = "no";
			if(  (chimeraResults[0].bsa >= minBS && chimeraResults[0].divr_qla_qrb >= divR)
			   ||
			   (chimeraResults[0].bsb >= minBS && chimeraResults[0].divr_qlb_qra >= divR) ) { chimeraFlag = "yes"; }
			
			
			if (chimeraFlag == "yes") {	
				if ((chimeraResults[0].bsa >= minBS) || (chimeraResults[0].bsb >= minBS)) {
					m->mothurOut(querySeq->getName() + "\tyes"); m->mothurOutEndLine();
					outAcc << querySeq->getName() << endl;
				}
			}
			
			printBlock(chimeraResults[0], out);
			out << endl;
		}else {  out << querySeq->getName() << "\tno" << endl;  }
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayer", "print");
		exit(1);
	}
}
//***************************************************************************************************************
int ChimeraSlayer::getChimeras(Sequence* query) {
	try {
		chimeraFlags = "no";
		
		//filter query
		spotMap = runFilter(query);
		
		querySeq = query;
		
		//referenceSeqs, numWanted, matchScore, misMatchPenalty, divR, minSimilarity
		maligner = new Maligner(templateSeqs, numWanted, match, misMatch, divR, minSim, minCov, searchMethod, databaseLeft, databaseRight);
		slayer = new Slayer(window, increment, minSim, divR, iters, minSNP);
		
		if (m->control_pressed) {  return 0;  }
		
		string chimeraFlag = maligner->getResults(query, decalc);
		if (m->control_pressed) {  return 0;  }
		vector<results> Results = maligner->getOutput();
				
		//found in testing realigning only made things worse
		if (realign) {
			ChimeraReAligner realigner(templateSeqs, match, misMatch);
			realigner.reAlign(query, Results);
		}

		if (chimeraFlag == "yes") {
			
			//get sequence that were given from maligner results
			vector<SeqDist> seqs;
			map<string, float> removeDups;
			map<string, float>::iterator itDup;
			map<string, string> parentNameSeq;
			map<string, string>::iterator itSeq;
			for (int j = 0; j < Results.size(); j++) {
				float dist = (Results[j].regionEnd - Results[j].regionStart + 1) * Results[j].queryToParentLocal;
				//only add if you are not a duplicate
				itDup = removeDups.find(Results[j].parent);
				if (itDup == removeDups.end()) { //this is not duplicate
					removeDups[Results[j].parent] = dist;
					parentNameSeq[Results[j].parent] = Results[j].parentAligned;
				}else if (dist > itDup->second) { //is this a stronger number for this parent
					removeDups[Results[j].parent] = dist;
					parentNameSeq[Results[j].parent] = Results[j].parentAligned;
				}
			}
			
			for (itDup = removeDups.begin(); itDup != removeDups.end(); itDup++) {
				//Sequence* seq = getSequence(itDup->first); //makes copy so you can filter and mask and not effect template
				itSeq = parentNameSeq.find(itDup->first);
//cout << itDup->first << itSeq->second << endl;
				Sequence* seq = new Sequence(itDup->first, itSeq->second);
				
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
			
			if (m->control_pressed) {  for (int k = 0; k < seqs.size(); k++) {  delete seqs[k].seq;   }  return 0;  }
			
			//send to slayer
			chimeraFlags = slayer->getResults(query, seqsForSlayer);
			if (m->control_pressed) {  return 0;  }
			chimeraResults = slayer->getOutput();
			
			//free memory
			for (int k = 0; k < seqs.size(); k++) {  delete seqs[k].seq;   }
		}
		
		delete maligner;
		delete slayer;
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayer", "getChimeras");
		exit(1);
	}
}
//***************************************************************************************************************
void ChimeraSlayer::printBlock(data_struct data, ostream& out){
	try {
	//out << "Name\tParentA\tParentB\tDivQLAQRB\tPerIDQLAQRB\tBootStrapA\tDivQLBQRA\tPerIDQLBQRA\tBootStrapB\tFlag\tLeftWindow\tRightWindow\n";
		
		out << querySeq->getName() << '\t';
		out << data.parentA.getName() << "\t" << data.parentB.getName()  << '\t';
		//out << "Left Window: " << spotMap[data.winLStart] << " " << spotMap[data.winLEnd] << endl;
		//out << "Right Window: " << spotMap[data.winRStart] << " " << spotMap[data.winREnd] << endl;
		
		out << data.divr_qla_qrb << '\t' << data.qla_qrb << '\t' << data.bsa << '\t';
		out << data.divr_qlb_qra << '\t' << data.qlb_qra << '\t' << data.bsb << '\t';
		
		out << "yes\t" << spotMap[data.winLStart] << "-" << spotMap[data.winLEnd] << '\t' << spotMap[data.winRStart] << "-" << spotMap[data.winREnd] << '\t';
		
		//out << "Similarity of parents: " << data.ab << endl;
		//out << "Similarity of query to parentA: " << data.qa << endl;
		//out << "Similarity of query to parentB: " << data.qb << endl;
		
		
		//out << "Per_id(QL,A): " << data.qla << endl;
		//out << "Per_id(QL,B): " << data.qlb << endl;
		//out << "Per_id(QR,A): " << data.qra << endl;
		//out << "Per_id(QR,B): " << data.qrb << endl;

		
		//out << "DeltaL: " << (data.qla - data.qlb) << endl;
		//out << "DeltaR: " << (data.qra - data.qrb) << endl;

	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayer", "printBlock");
		exit(1);
	}
}
//***************************************************************************************************************

