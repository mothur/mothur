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
#include "blastdb.hpp"

//***************************************************************************************************************
ChimeraSlayer::ChimeraSlayer(string file, string temp, bool trim, string mode, int k, int ms, int mms, int win, float div, 
int minsim, int mincov, int minbs, int minsnp, int par, int it, int inc, int numw, bool r, string blas, int tid) : Chimera()  {  	
	try {
		fastafile = file;
		templateFileName = temp; templateSeqs = readSeqs(temp);
		searchMethod = mode;
		kmerSize = k;
		match = ms;
		misMatch = mms;
		window = win;
		divR = div;
		minSim = minsim;
		minCov = mincov;
		minBS = minbs;
		minSNP = minsnp;
		parents = par;
		iters = it;
		increment = inc;
		numWanted = numw;
		realign = r; 
		trimChimera = trim;
		numNoParents = 0;
		blastlocation = blas;
		threadID = tid;
	
		doPrep();
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayer", "ChimeraSlayer");
		exit(1);
	}
}
//***************************************************************************************************************
//template=self
ChimeraSlayer::ChimeraSlayer(string file, string temp, bool trim, map<string, int>& prior, string mode, int k, int ms, int mms, int win, float div, 
							 int minsim, int mincov, int minbs, int minsnp, int par, int it, int inc, int numw, bool r, string blas, int tid, bool bg) : Chimera()  {  	
	try {
		byGroup = bg;
		fastafile = file; templateSeqs = readSeqs(fastafile);
		templateFileName = temp; 
		searchMethod = mode;
		kmerSize = k;
		match = ms;
		misMatch = mms;
		window = win;
		divR = div;
		minSim = minsim;
		minCov = mincov;
		minBS = minbs;
		minSNP = minsnp;
		parents = par;
		iters = it;
		increment = inc;
		numWanted = numw;
		realign = r; 
		trimChimera = trim;
		priority = prior;
		numNoParents = 0;
		blastlocation = blas;
		threadID = tid;
		
		
		createFilter(templateSeqs, 0.0); //just removed columns where all seqs have a gap
		
		if (searchMethod == "distance") { 
			//createFilter(templateSeqs, 0.0); //just removed columns where all seqs have a gap
			
			//run filter on template copying templateSeqs into filteredTemplateSeqs
			for (int i = 0; i < templateSeqs.size(); i++) {  
				if (m->control_pressed) {  break; }
				
				Sequence* newSeq = new Sequence(templateSeqs[i]->getName(), templateSeqs[i]->getAligned());
				runFilter(newSeq);  
				filteredTemplateSeqs.push_back(newSeq);
			}
		}
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayer", "ChimeraSlayer");
		exit(1);
	}
}
//***************************************************************************************************************
//template=self
ChimeraSlayer::ChimeraSlayer(string file, string temp, bool trim, map<string, int>& prior, string mode, int k, int ms, int mms, int win, float div, 
							 int minsim, int mincov, int minbs, int minsnp, int par, int it, int inc, int numw, bool r, string blas, int tid) : Chimera()  {  	
	try {
		fastafile = file; templateSeqs = readSeqs(fastafile);
		templateFileName = temp; 
		searchMethod = mode;
		kmerSize = k;
		match = ms;
		misMatch = mms;
		window = win;
		divR = div;
		minSim = minsim;
		minCov = mincov;
		minBS = minbs;
		minSNP = minsnp;
		parents = par;
		iters = it;
		increment = inc;
		numWanted = numw;
		realign = r; 
		trimChimera = trim;
		priority = prior;
		numNoParents = 0;
		blastlocation = blas;
		threadID = tid;
		
		
		createFilter(templateSeqs, 0.0); //just removed columns where all seqs have a gap
		
		if (searchMethod == "distance") { 
			//createFilter(templateSeqs, 0.0); //just removed columns where all seqs have a gap
			
			//run filter on template copying templateSeqs into filteredTemplateSeqs
			for (int i = 0; i < templateSeqs.size(); i++) {  
				if (m->control_pressed) {  break; }
				
				Sequence* newSeq = new Sequence(templateSeqs[i]->getName(), templateSeqs[i]->getAligned());
				runFilter(newSeq);  
				filteredTemplateSeqs.push_back(newSeq);
			}
		}
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayer", "ChimeraSlayer");
		exit(1);
	}
}
//***************************************************************************************************************
int ChimeraSlayer::doPrep() {
	try {
		if (searchMethod == "distance") { 
			//read in all query seqs
			vector<Sequence*> tempQuerySeqs = readSeqs(fastafile);
		
			vector<Sequence*> temp = templateSeqs;
			for (int i = 0; i < tempQuerySeqs.size(); i++) {  temp.push_back(tempQuerySeqs[i]);  }
		
			createFilter(temp, 0.0); //just removed columns where all seqs have a gap
		
			for (int i = 0; i < tempQuerySeqs.size(); i++) { delete tempQuerySeqs[i];  }
		
			if (m->control_pressed) {  return 0; } 
		
			//run filter on template copying templateSeqs into filteredTemplateSeqs
			for (int i = 0; i < templateSeqs.size(); i++) {  
				if (m->control_pressed) {  return 0; }
				
				Sequence* newSeq = new Sequence(templateSeqs[i]->getName(), templateSeqs[i]->getAligned());
				runFilter(newSeq);  
				filteredTemplateSeqs.push_back(newSeq);
			}
		}
		string 	kmerDBNameLeft;
		string 	kmerDBNameRight;
	
		//generate the kmerdb to pass to maligner
		if (searchMethod == "kmer") { 
			string templatePath = m->hasPath(templateFileName);
			string rightTemplateFileName = templatePath + "right." + m->getRootName(m->getSimpleName(templateFileName));
			databaseRight = new KmerDB(rightTemplateFileName, kmerSize);
				
			string leftTemplateFileName = templatePath + "left." + m->getRootName(m->getSimpleName(templateFileName));
			databaseLeft = new KmerDB(leftTemplateFileName, kmerSize);	
		
			//leftside
			kmerDBNameLeft = leftTemplateFileName.substr(0,leftTemplateFileName.find_last_of(".")+1) + char('0'+ kmerSize) + "mer";
			ifstream kmerFileTestLeft(kmerDBNameLeft.c_str());
			bool needToGenerateLeft = true;
			
			if(kmerFileTestLeft){	
				bool GoodFile = m->checkReleaseVersion(kmerFileTestLeft, m->getVersion());
				if (GoodFile) {  needToGenerateLeft = false;	}
			}
			
			if(needToGenerateLeft){	
			
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
			kmerDBNameRight = rightTemplateFileName.substr(0,rightTemplateFileName.find_last_of(".")+1) + char('0'+ kmerSize) + "mer";
			ifstream kmerFileTestRight(kmerDBNameRight.c_str());
			bool needToGenerateRight = true;
			
			if(kmerFileTestRight){	
				bool GoodFile = m->checkReleaseVersion(kmerFileTestRight, m->getVersion());
				if (GoodFile) {  needToGenerateRight = false;	}
			}
			
			if(needToGenerateRight){	
			
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
    
		}else if (searchMethod == "blast") {
		
			//generate blastdb
			databaseLeft = new BlastDB(m->getRootName(m->getSimpleName(fastafile)), -1.0, -1.0, 1, -3, blastlocation, threadID);
			
			if (m->control_pressed) { return 0; }

			for (int i = 0; i < templateSeqs.size(); i++) { 	databaseLeft->addSequence(*templateSeqs[i]);	}
			databaseLeft->generateDB();
			databaseLeft->setNumSeqs(templateSeqs.size());
		}
		
		return 0;

	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayer", "doprep");
		exit(1);
	}
}
//***************************************************************************************************************
vector<Sequence*> ChimeraSlayer::getTemplate(Sequence q, vector<Sequence*>& userTemplateFiltered) {
	try {
		
		//when template=self, the query file is sorted from most abundance to least abundant
		//userTemplate grows as the query file is processed by adding sequences that are not chimeric and more abundant
		vector<Sequence*> userTemplate;
		
		int myAbund = priority[q.getName()];
		
		for (int i = 0; i < templateSeqs.size(); i++) {
			
			if (m->control_pressed) { return userTemplate; } 
			
			//have I reached a sequence with the same abundance as myself?
			if (!(priority[templateSeqs[i]->getName()] > myAbund)) { break; }
			
			//if its am not chimeric add it
			if (chimericSeqs.count(templateSeqs[i]->getName()) == 0) { 
				userTemplate.push_back(templateSeqs[i]); 
				if (searchMethod == "distance") { userTemplateFiltered.push_back(filteredTemplateSeqs[i]); }
			}
		}
		
		//avoids nuisance error from formatdb for making blank blast database
		if (userTemplate.size() == 0) {
			return userTemplate;
		}
		
		string 	kmerDBNameLeft;
		string 	kmerDBNameRight;
		
		//generate the kmerdb to pass to maligner
		if (searchMethod == "kmer") { 
			string templatePath = m->hasPath(templateFileName);
			string rightTemplateFileName = templatePath + "right." + m->getRootName(m->getSimpleName(templateFileName));
			databaseRight = new KmerDB(rightTemplateFileName, kmerSize);
			
			string leftTemplateFileName = templatePath + "left." + m->getRootName(m->getSimpleName(templateFileName));
			databaseLeft = new KmerDB(leftTemplateFileName, kmerSize);	
			
			
			for (int i = 0; i < userTemplate.size(); i++) {
				
				if (m->control_pressed) { return userTemplate; } 
				
				string leftFrag = userTemplate[i]->getUnaligned();
				leftFrag = leftFrag.substr(0, int(leftFrag.length() * 0.33));
				
				Sequence leftTemp(userTemplate[i]->getName(), leftFrag);
				databaseLeft->addSequence(leftTemp);	
			}
			databaseLeft->generateDB();
			databaseLeft->setNumSeqs(userTemplate.size());
				
			for (int i = 0; i < userTemplate.size(); i++) {
				if (m->control_pressed) { return userTemplate; }  
					
				string rightFrag = userTemplate[i]->getUnaligned();
				rightFrag = rightFrag.substr(int(rightFrag.length() * 0.66));
					
				Sequence rightTemp(userTemplate[i]->getName(), rightFrag);
				databaseRight->addSequence(rightTemp);	
			}
			databaseRight->generateDB();
			databaseRight->setNumSeqs(userTemplate.size());
	
		}else if (searchMethod == "blast") {
			
			//generate blastdb
			databaseLeft = new BlastDB(m->getRootName(m->getSimpleName(templateFileName)), -1.0, -1.0, 1, -3, blastlocation, threadID);
			
			if (m->control_pressed) { return userTemplate; }

			for (int i = 0; i < userTemplate.size(); i++) { if (m->control_pressed) { return userTemplate; }   databaseLeft->addSequence(*userTemplate[i]);	}
			databaseLeft->generateDB();
			databaseLeft->setNumSeqs(userTemplate.size());
		}
		
		return userTemplate;
		
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayer", "getTemplate");
		exit(1);
	}
}

//***************************************************************************************************************
ChimeraSlayer::~ChimeraSlayer() { 	
	if (templateFileName != "self") {
		if (searchMethod == "kmer") {  delete databaseRight;  delete databaseLeft;  }	
		else if (searchMethod == "blast") {  delete databaseLeft; }
	}
}
//***************************************************************************************************************
void ChimeraSlayer::printHeader(ostream& out) {
	m->mothurOutEndLine();
	m->mothurOut("Only reporting sequence supported by " + toString(minBS) + "% of bootstrapped results.");
	m->mothurOutEndLine();
	
	out << "Name\tLeftParent\tRightParent\tDivQLAQRB\tPerIDQLAQRB\tBootStrapA\tDivQLBQRA\tPerIDQLBQRA\tBootStrapB\tFlag\tLeftWindow\tRightWindow\n";
}
//***************************************************************************************************************
Sequence ChimeraSlayer::print(ostream& out, ostream& outAcc) {
	try {
		Sequence trim;
		if (trimChimera) { trim.setName(trimQuery.getName()); trim.setAligned(trimQuery.getAligned()); }
		
		if (chimeraFlags == "yes") {
			string chimeraFlag = "no";
			if(  (chimeraResults[0].bsa >= minBS && chimeraResults[0].divr_qla_qrb >= divR)
			   ||
			   (chimeraResults[0].bsb >= minBS && chimeraResults[0].divr_qlb_qra >= divR) ) { chimeraFlag = "yes"; }
			
			
			if (chimeraFlag == "yes") {	
				if ((chimeraResults[0].bsa >= minBS) || (chimeraResults[0].bsb >= minBS)) {
					m->mothurOut(querySeq.getName() + "\tyes"); m->mothurOutEndLine();
					outAcc << querySeq.getName() << endl;
					
					if (templateFileName == "self") {  chimericSeqs.insert(querySeq.getName()); }
					
					if (trimChimera) {  
						int lengthLeft = chimeraResults[0].winLEnd - chimeraResults[0].winLStart;
						int lengthRight = chimeraResults[0].winREnd - chimeraResults[0].winRStart;
						
						string newAligned = trim.getAligned();

						if (lengthLeft > lengthRight) { //trim right
							for (int i = (chimeraResults[0].winRStart-1); i < newAligned.length(); i++) { newAligned[i] = '.'; }
						}else { //trim left
							for (int i = 0; i < chimeraResults[0].winLEnd; i++) { newAligned[i] = '.'; }
						}
						trim.setAligned(newAligned);
					}
				}
			}
			
			printBlock(chimeraResults[0], chimeraFlag, out);
		}else {  
			out << querySeq.getName() << "\tno" << endl; 
		}
		
		return trim;
		
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayer", "print");
		exit(1);
	}
}
//***************************************************************************************************************
Sequence ChimeraSlayer::print(ostream& out, ostream& outAcc, data_results leftPiece, data_results rightPiece) {
	try {
		Sequence trim;
				
		if (trimChimera) { 
			string aligned = leftPiece.trimQuery.getAligned() + rightPiece.trimQuery.getAligned();
			trim.setName(leftPiece.trimQuery.getName()); trim.setAligned(aligned); 
		}
		
		if ((leftPiece.flag == "yes") || (rightPiece.flag == "yes")) {
			
			string chimeraFlag = "no";
			if (leftPiece.flag == "yes") {
				
				if(  (leftPiece.results[0].bsa >= minBS && leftPiece.results[0].divr_qla_qrb >= divR)
					||
					(leftPiece.results[0].bsb >= minBS && leftPiece.results[0].divr_qlb_qra >= divR) ) { chimeraFlag = "yes"; }
			}
			
			if (rightPiece.flag == "yes") {
				if ( (rightPiece.results[0].bsa >= minBS && rightPiece.results[0].divr_qla_qrb >= divR)
				 ||
				 (rightPiece.results[0].bsb >= minBS && rightPiece.results[0].divr_qlb_qra >= divR) ) { chimeraFlag = "yes"; }
			}
			
			bool rightChimeric = false;
			bool leftChimeric = false;
			
			if (chimeraFlag == "yes") {	
				//which peice is chimeric or are both
				if (rightPiece.flag == "yes") { if ((rightPiece.results[0].bsa >= minBS) || (rightPiece.results[0].bsb >= minBS)) { rightChimeric = true; } }
				if (leftPiece.flag == "yes") { if ((leftPiece.results[0].bsa >= minBS) || (leftPiece.results[0].bsb >= minBS))	{ leftChimeric = true;	} }
				
				if (rightChimeric || leftChimeric) {
					m->mothurOut(querySeq.getName() + "\tyes"); m->mothurOutEndLine();
					outAcc << querySeq.getName() << endl;
					
					if (templateFileName == "self") {  chimericSeqs.insert(querySeq.getName()); }
					
					if (trimChimera) {  
						string newAligned = trim.getAligned();
												
						//right side is fine so keep that
						if ((leftChimeric) && (!rightChimeric)) {
							for (int i = 0; i < leftPiece.results[0].winREnd; i++) { newAligned[i] = '.'; } 
						}else if ((!leftChimeric) && (rightChimeric)) { //leftside is fine so keep that
							for (int i = (rightPiece.results[0].winLStart-1); i < newAligned.length(); i++) { newAligned[i] = '.'; }
						}else { //both sides are chimeric, keep longest piece
							
							int lengthLeftLeft = leftPiece.results[0].winLEnd - leftPiece.results[0].winLStart;
							int lengthLeftRight = leftPiece.results[0].winREnd - leftPiece.results[0].winRStart;
							
							int longest = 1; // leftleft = 1, leftright = 2, rightleft = 3 rightright = 4
							int length = lengthLeftLeft;
							if (lengthLeftLeft < lengthLeftRight) { longest = 2;  length = lengthLeftRight; }
							
							int lengthRightLeft = rightPiece.results[0].winLEnd - rightPiece.results[0].winLStart;
							int lengthRightRight = rightPiece.results[0].winREnd - rightPiece.results[0].winRStart;
							
							if (lengthRightLeft > length) { longest = 3; length = lengthRightLeft;  }
							if (lengthRightRight > length) { longest = 4; }
							
							if (longest == 1) { //leftleft
								for (int i = (leftPiece.results[0].winRStart-1); i < newAligned.length(); i++) { newAligned[i] = '.'; }
							}else if (longest == 2) { //leftright
								//get rid of leftleft
								for (int i = (leftPiece.results[0].winLStart-1); i < (leftPiece.results[0].winLEnd-1); i++) { newAligned[i] = '.'; }
								//get rid of right
								for (int i = (rightPiece.results[0].winLStart-1); i < newAligned.length(); i++) { newAligned[i] = '.'; }
							}else if (longest == 3) { //rightleft
								//get rid of left
								for (int i = 0; i < leftPiece.results[0].winREnd; i++) { newAligned[i] = '.'; } 
								//get rid of rightright
								for (int i = (rightPiece.results[0].winRStart-1); i < newAligned.length(); i++) { newAligned[i] = '.'; }
							}else { //rightright
								//get rid of left
								for (int i = 0; i < leftPiece.results[0].winREnd; i++) { newAligned[i] = '.'; } 
								//get rid of rightleft
								for (int i = (rightPiece.results[0].winLStart-1); i < (rightPiece.results[0].winLEnd-1); i++) { newAligned[i] = '.'; }
							}
						}
							
						trim.setAligned(newAligned);
					}
					
				}
			}
			
			printBlock(leftPiece, rightPiece, leftChimeric, rightChimeric, chimeraFlag, out);
		}else {  
			out << querySeq.getName() << "\tno" << endl;  
		}
		
		return trim;
		
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayer", "print");
		exit(1);
	}
}

//***************************************************************************************************************
int ChimeraSlayer::getChimeras(Sequence* query) {
	try {
		
		trimQuery.setName(query->getName()); trimQuery.setAligned(query->getAligned());
		printResults.trimQuery = trimQuery; 
		
		chimeraFlags = "no";
		printResults.flag = "no";
		
		querySeq = *query;
		
		//you must create a template
		vector<Sequence*> thisTemplate;
		vector<Sequence*> thisFilteredTemplate;
		if (templateFileName != "self") { thisTemplate = templateSeqs; thisFilteredTemplate = filteredTemplateSeqs; }
		else {  thisTemplate = getTemplate(*query, thisFilteredTemplate);  } //fills this template and creates the databases
		
		if (m->control_pressed) {  return 0;  }
		if (thisTemplate.size() == 0) {  return 0; } //not chimeric
		
		//moved this out of maligner - 4/29/11
		vector<Sequence> refSeqs = getRefSeqs(*query, thisTemplate, thisFilteredTemplate);
		
		Maligner maligner(refSeqs, match, misMatch, divR, minSim, minCov); 
		Slayer slayer(window, increment, minSim, divR, iters, minSNP, minBS);
		
		if (templateFileName == "self") {
			if (searchMethod == "kmer") {  delete databaseRight;  delete databaseLeft;  }	
			else if (searchMethod == "blast") {  delete databaseLeft; }
		}
	
		if (m->control_pressed) {  return 0;  }

		string chimeraFlag = maligner.getResults(*query, decalc);

		if (m->control_pressed) {  return 0;  }
		
		vector<results> Results = maligner.getOutput();
		
		//for (int i = 0; i < refSeqs.size(); i++) {  delete refSeqs[i];	}
		
		if (chimeraFlag == "yes") {
			
			if (realign) {
				vector<string> parents;
				for (int i = 0; i < Results.size(); i++) {
					parents.push_back(Results[i].parentAligned);
				}
				
				ChimeraReAligner realigner;		
				realigner.reAlign(query, parents);

			}
			
//			cout << query->getAligned() << endl;
			//get sequence that were given from maligner results
			vector<SeqCompare> seqs;
			map<string, float> removeDups;
			map<string, float>::iterator itDup;
			map<string, string> parentNameSeq;
			map<string, string>::iterator itSeq;
			for (int j = 0; j < Results.size(); j++) {

				float dist = (Results[j].regionEnd - Results[j].regionStart + 1) * Results[j].queryToParentLocal;
				//only add if you are not a duplicate
//				cout << Results[j].parent << '\t' << Results[j].regionEnd << '\t' << Results[j].regionStart << '\t' << Results[j].regionEnd - Results[j].regionStart +1 << '\t' << Results[j].queryToParentLocal << '\t' << dist << endl;
				
				
				if(Results[j].queryToParentLocal >= 90){	//local match has to be over 90% similarity
				
					itDup = removeDups.find(Results[j].parent);
					if (itDup == removeDups.end()) { //this is not duplicate
						removeDups[Results[j].parent] = dist;
						parentNameSeq[Results[j].parent] = Results[j].parentAligned;
					}else if (dist > itDup->second) { //is this a stronger number for this parent
						removeDups[Results[j].parent] = dist;
						parentNameSeq[Results[j].parent] = Results[j].parentAligned;
					}
				
				}
				
			}
			
			for (itDup = removeDups.begin(); itDup != removeDups.end(); itDup++) {
				itSeq = parentNameSeq.find(itDup->first);
				Sequence seq(itDup->first, itSeq->second);
				
				SeqCompare member;
				member.seq = seq;
				member.dist = itDup->second;
				seqs.push_back(member);
			}
			
			//limit number of parents to explore - default 3
			if (Results.size() > parents) {
				//sort by distance
				sort(seqs.begin(), seqs.end(), compareSeqCompare);
				//prioritize larger more similiar sequence fragments
				reverse(seqs.begin(), seqs.end());
				
				//for (int k = seqs.size()-1; k > (parents-1); k--)  {  
				//	delete seqs[k].seq;
					//seqs.pop_back();	
				//}
			}
		
			//put seqs into vector to send to slayer
			
//			cout << query->getAligned() << endl;
			vector<Sequence> seqsForSlayer;
			for (int k = 0; k < seqs.size(); k++) {  
//				cout << seqs[k].seq->getAligned() << endl;
				seqsForSlayer.push_back(seqs[k].seq);	
//				cout << seqs[k].seq->getName() << endl;
			}
			
			if (m->control_pressed) {  return 0;  }

			//send to slayer
			chimeraFlags = slayer.getResults(*query, seqsForSlayer);
			if (m->control_pressed) {  return 0;  }
			chimeraResults = slayer.getOutput();
			
			printResults.flag = chimeraFlags;
			printResults.results = chimeraResults;
			
			//free memory
			//for (int k = 0; k < seqs.size(); k++) {  delete seqs[k].seq;   }
		}
		//cout << endl << endl;
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayer", "getChimeras");
		exit(1);
	}
}
//***************************************************************************************************************
void ChimeraSlayer::printBlock(data_struct data, string flag, ostream& out){
	try {
		out << querySeq.getName();
		out << '\t' << data.parentA.getName() << "\t" << data.parentB.getName();
	
		out << '\t' <<  data.divr_qla_qrb << '\t' << data.qla_qrb << '\t' << data.bsa;
		out << '\t' <<  data.divr_qlb_qra << '\t' << data.qlb_qra << '\t' << data.bsb ;
		
		out << '\t' <<  flag << '\t' << data.winLStart << "-" << data.winLEnd << '\t' << data.winRStart << "-" << data.winREnd << '\n';
		
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayer", "printBlock");
		exit(1);
	}
}
//***************************************************************************************************************
void ChimeraSlayer::printBlock(data_results leftdata, data_results rightdata, bool leftChimeric, bool rightChimeric, string flag, ostream& out){
	try {
		
		if ((leftChimeric) && (!rightChimeric)) { //print left
			out << querySeq.getName();
			out << '\t' << leftdata.results[0].parentA.getName() << "\t" << leftdata.results[0].parentB.getName();
			out << '\t' << leftdata.results[0].divr_qla_qrb << '\t' << leftdata.results[0].qla_qrb << '\t' << leftdata.results[0].bsa;
			out << '\t' << leftdata.results[0].divr_qlb_qra << '\t' << leftdata.results[0].qlb_qra << '\t' << leftdata.results[0].bsb;
            out << '\t' << flag << '\t' << leftdata.results[0].winLStart << "-" << leftdata.results[0].winLEnd << '\t' << leftdata.results[0].winRStart << "-" << leftdata.results[0].winREnd << endl;
		
		}else if ((!leftChimeric) && (rightChimeric)) {  //print right
			out << querySeq.getName();
			out << '\t' << rightdata.results[0].parentA.getName() << "\t" << rightdata.results[0].parentB.getName();
			out << '\t' << rightdata.results[0].divr_qla_qrb << '\t' << rightdata.results[0].qla_qrb << '\t' << rightdata.results[0].bsa;
			out << '\t' << rightdata.results[0].divr_qlb_qra << '\t' << rightdata.results[0].qlb_qra << '\t' << rightdata.results[0].bsb;
			out << '\t' << flag << '\t' << rightdata.results[0].winLStart << "-" << rightdata.results[0].winLEnd << '\t' << rightdata.results[0].winRStart << "-" << rightdata.results[0].winREnd << endl;
			
		}else  { //print both results
			if (leftdata.flag == "yes") {
				out << querySeq.getName() + "_LEFT";
				out << '\t' << leftdata.results[0].parentA.getName() << "\t" << leftdata.results[0].parentB.getName();
				out << '\t' << leftdata.results[0].divr_qla_qrb << '\t' << leftdata.results[0].qla_qrb << '\t' << leftdata.results[0].bsa;
				out << '\t' << leftdata.results[0].divr_qlb_qra << '\t' << leftdata.results[0].qlb_qra << '\t' << leftdata.results[0].bsb;
                out << '\t' << flag << '\t' << leftdata.results[0].winLStart << "-" << leftdata.results[0].winLEnd << '\t' << leftdata.results[0].winRStart << "-" << leftdata.results[0].winREnd << endl;
			}
			
			if (rightdata.flag == "yes") {
				out << querySeq.getName() + "_RIGHT";
				out << '\t' << rightdata.results[0].parentA.getName() << "\t" << rightdata.results[0].parentB.getName();
				out << '\t' << rightdata.results[0].divr_qla_qrb << '\t' << rightdata.results[0].qla_qrb << '\t' << rightdata.results[0].bsa;
				out << '\t' << rightdata.results[0].divr_qlb_qra << '\t' << rightdata.results[0].qlb_qra << '\t' << rightdata.results[0].bsb;
				out << '\t' << flag << '\t' << rightdata.results[0].winLStart << "-" << rightdata.results[0].winLEnd << '\t' << rightdata.results[0].winRStart << "-" << rightdata.results[0].winREnd << '\n';
		
			}
		}
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayer", "printBlock");
		exit(1);
	}
}
//***************************************************************************************************************
string ChimeraSlayer::getBlock(data_results leftdata, data_results rightdata, bool leftChimeric, bool rightChimeric, string flag){
	try {
		
		string out = "";
		
		if ((leftChimeric) && (!rightChimeric)) { //get left
			out += querySeq.getName();
			out += "\t" + leftdata.results[0].parentA.getName() + "\t" + leftdata.results[0].parentB.getName();
			out += "\t" + toString(leftdata.results[0].divr_qla_qrb) + "\t" + toString(leftdata.results[0].qla_qrb) + "\t" + toString(leftdata.results[0].bsa);
			out += "\t" + toString(leftdata.results[0].divr_qlb_qra) + "\t" + toString(leftdata.results[0].qlb_qra) + "\t" + toString(leftdata.results[0].bsb);
			out += "\t" + flag + "\t" + toString(leftdata.results[0].winLStart) + "-" + toString(leftdata.results[0].winLEnd) + "\t" + toString(leftdata.results[0].winRStart) + "-" + toString(leftdata.results[0].winREnd) + "\n";
			
		}else if ((!leftChimeric) && (rightChimeric)) {  //print right
			out += querySeq.getName();
			out += "\t" + rightdata.results[0].parentA.getName() + "\t" + rightdata.results[0].parentB.getName();
			out += "\t" + toString(rightdata.results[0].divr_qla_qrb) + "\t" + toString(rightdata.results[0].qla_qrb) + "\t" + toString(rightdata.results[0].bsa);
			out += "\t" + toString(rightdata.results[0].divr_qlb_qra) + "\t" + toString(rightdata.results[0].qlb_qra) + "\t" + toString(rightdata.results[0].bsb);
			out += "\t" + flag + "\t" + toString(rightdata.results[0].winLStart) + "-" + toString(rightdata.results[0].winLEnd) + "\t" + toString(rightdata.results[0].winRStart) + "-" + toString(rightdata.results[0].winREnd) + "\n";
			
		}else  { //print both results
			if (leftdata.flag == "yes") {
				out += querySeq.getName() + "_LEFT";
				out += "\t" + leftdata.results[0].parentA.getName() + "\t" + leftdata.results[0].parentB.getName();
				out += "\t" + toString(leftdata.results[0].divr_qla_qrb) + "\t" + toString(leftdata.results[0].qla_qrb) + "\t" + toString(leftdata.results[0].bsa);
				out += "\t" + toString(leftdata.results[0].divr_qlb_qra) + "\t" + toString(leftdata.results[0].qlb_qra) + "\t" + toString(leftdata.results[0].bsb);
				out += "\t" + flag + "\t" + toString(leftdata.results[0].winLStart) + "-" + toString(leftdata.results[0].winLEnd) + "\t" + toString(leftdata.results[0].winRStart) + "-" + toString(leftdata.results[0].winREnd) + "\n";
			}
			
			if (rightdata.flag == "yes") {
				out +=  querySeq.getName() + "_RIGHT";
				out += "\t" + rightdata.results[0].parentA.getName() + "\t" + rightdata.results[0].parentB.getName();
				out += "\t" + toString(rightdata.results[0].divr_qla_qrb) + "\t" + toString(rightdata.results[0].qla_qrb) + "\t" + toString(rightdata.results[0].bsa);
				out += "\t" + toString(rightdata.results[0].divr_qlb_qra) + "\t" + toString(rightdata.results[0].qlb_qra) + "\t" + toString(rightdata.results[0].bsb);
				out += "\t" + flag + "\t" + toString(rightdata.results[0].winLStart) + "-" + toString(rightdata.results[0].winLEnd) + "\t" + toString(rightdata.results[0].winRStart) + "-" + toString(rightdata.results[0].winREnd) + "\n";
			}
		}
		
		return out;
		
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayer", "getBlock");
		exit(1);
	}
}
//***************************************************************************************************************
string ChimeraSlayer::getBlock(data_struct data, string flag){
	try {
		
		string outputString = "";
		outputString += querySeq.getName();
		outputString +=  "\t" + data.parentA.getName() + "\t" + data.parentB.getName();
		outputString +=  "\t" + toString(data.divr_qla_qrb) + "\t" + toString(data.qla_qrb) + "\t" + toString(data.bsa);
		outputString +=  "\t" + toString(data.divr_qlb_qra) + "\t" + toString(data.qlb_qra) + "\t" + toString(data.bsb);
		outputString +=  "\t" + flag + "\t" + toString(data.winLStart) + "-" + toString(data.winLEnd) + "\t" + toString(data.winRStart) + "-" + toString(data.winREnd) + "\n";
		
		return outputString;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayer", "getBlock");
		exit(1);
	}
}
//***************************************************************************************************************
vector<Sequence> ChimeraSlayer::getRefSeqs(Sequence q, vector<Sequence*>& thisTemplate, vector<Sequence*>& thisFilteredTemplate){
	try {
		
		vector<Sequence> refSeqs;
		
		if (searchMethod == "distance") {
			//find closest seqs to query in template - returns copies of seqs so trim does not destroy - remember to deallocate
			Sequence* newSeq = new Sequence(q.getName(), q.getAligned());
			runFilter(newSeq);
			refSeqs = decalc.findClosest(*newSeq, thisTemplate, thisFilteredTemplate, numWanted, minSim);
			delete newSeq;
		}else if (searchMethod == "blast")  {
			refSeqs = getBlastSeqs(q, thisTemplate, numWanted); //fills indexes
		}else if (searchMethod == "kmer") {
			refSeqs = getKmerSeqs(q, thisTemplate, numWanted); //fills indexes
		}else { m->mothurOut("not valid search."); exit(1);  } //should never get here
		
		return refSeqs;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayer", "getRefSeqs");
		exit(1);
	}
}
//***************************************************************************************************************/
vector<Sequence> ChimeraSlayer::getBlastSeqs(Sequence q, vector<Sequence*>& db, int num) {
	try {	
		
		vector<Sequence> refResults;
		
		//get parts of query
		string queryUnAligned = q.getUnaligned();
		string leftQuery = queryUnAligned.substr(0, int(queryUnAligned.length() * 0.33)); //first 1/3 of the sequence
		string rightQuery = queryUnAligned.substr(int(queryUnAligned.length() * 0.66)); //last 1/3 of the sequence
//cout << "whole length = " << queryUnAligned.length() << '\t' << "left length = " << leftQuery.length() << '\t' << "right length = "<< rightQuery.length() << endl;	
		Sequence* queryLeft = new Sequence(q.getName(), leftQuery);
		Sequence* queryRight = new Sequence(q.getName(), rightQuery);
		
		vector<int> tempIndexesLeft = databaseLeft->findClosestMegaBlast(queryLeft, num+1, minSim);
		vector<int> tempIndexesRight = databaseLeft->findClosestMegaBlast(queryRight, num+1, minSim);
				
		
		//cout << q->getName() << '\t' << leftQuery << '\t' << "leftMatches = " << tempIndexesLeft.size() << '\t' << rightQuery	<< " rightMatches = " << tempIndexesRight.size() << endl;
//		vector<int> smaller;
//		vector<int> larger;
//		
//		if (tempIndexesRight.size() < tempIndexesLeft.size()) { smaller = tempIndexesRight;  larger = tempIndexesLeft;  }
//		else { smaller = tempIndexesLeft;  larger = tempIndexesRight;  } 
		
		//merge results		
		map<int, int> seen;
		map<int, int>::iterator it;
		vector<int> mergedResults;
		
		int index = 0;
//		for (int i = 0; i < smaller.size(); i++) {
		while(index < tempIndexesLeft.size() && index < tempIndexesRight.size()){
			
			if (m->control_pressed) { delete queryRight; delete queryLeft; return refResults; }
	
			//add left if you havent already
			it = seen.find(tempIndexesLeft[index]);
			if (it == seen.end()) {  
				mergedResults.push_back(tempIndexesLeft[index]);
				seen[tempIndexesLeft[index]] = tempIndexesLeft[index];
			}
			
			//add right if you havent already
			it = seen.find(tempIndexesRight[index]);
			if (it == seen.end()) {  
				mergedResults.push_back(tempIndexesRight[index]);
				seen[tempIndexesRight[index]] = tempIndexesRight[index];
			}
			index++;
		}

		
		for (int i = index; i < tempIndexesLeft.size(); i++) {
			if (m->control_pressed) { delete queryRight; delete queryLeft; return refResults; }
			
			//add right if you havent already
			it = seen.find(tempIndexesLeft[i]);
			if (it == seen.end()) {  
				mergedResults.push_back(tempIndexesLeft[i]);
				seen[tempIndexesLeft[i]] = tempIndexesLeft[i];
			}
		}

		for (int i = index; i < tempIndexesRight.size(); i++) {
			if (m->control_pressed) { delete queryRight; delete queryLeft; return refResults; }
			
			//add right if you havent already
			it = seen.find(tempIndexesRight[i]);
			if (it == seen.end()) {  
				mergedResults.push_back(tempIndexesRight[i]);
				seen[tempIndexesRight[i]] = tempIndexesRight[i];
			}
		}
		//string qname = q->getName().substr(0, q->getName().find_last_of('_'));	
		//cout << qname << endl;	
		
		if (mergedResults.size() == 0) { numNoParents++; }
		
		for (int i = 0; i < mergedResults.size(); i++) {
			//cout << q->getName() << mergedResults[i]  << '\t' << db[mergedResults[i]]->getName() << endl;	
			if (db[mergedResults[i]]->getName() != q.getName()) { 
				Sequence temp(db[mergedResults[i]]->getName(), db[mergedResults[i]]->getAligned());
				refResults.push_back(temp);
			}
		}
		//cout << endl << endl;

		delete queryRight;
		delete queryLeft;
		
		return refResults;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayer", "getBlastSeqs");
		exit(1);
	}
}
//***************************************************************************************************************
vector<Sequence> ChimeraSlayer::getKmerSeqs(Sequence q, vector<Sequence*>& db, int num) {
	try {	
		vector<Sequence> refResults;
		
		//get parts of query
		string queryUnAligned = q.getUnaligned();
		string leftQuery = queryUnAligned.substr(0, int(queryUnAligned.length() * 0.33)); //first 1/3 of the sequence
		string rightQuery = queryUnAligned.substr(int(queryUnAligned.length() * 0.66)); //last 1/3 of the sequence
		
		Sequence* queryLeft = new Sequence(q.getName(), leftQuery);
		Sequence* queryRight = new Sequence(q.getName(), rightQuery);
		
		vector<int> tempIndexesLeft = databaseLeft->findClosestSequences(queryLeft, num);
		vector<int> tempIndexesRight = databaseRight->findClosestSequences(queryRight, num);
		
		//merge results		
		map<int, int> seen;
		map<int, int>::iterator it;
		vector<int> mergedResults;
		
		int index = 0;
		//		for (int i = 0; i < smaller.size(); i++) {
		while(index < tempIndexesLeft.size() && index < tempIndexesRight.size()){
			
			if (m->control_pressed) { delete queryRight; delete queryLeft; return refResults; }
			
			//add left if you havent already
			it = seen.find(tempIndexesLeft[index]);
			if (it == seen.end()) {  
				mergedResults.push_back(tempIndexesLeft[index]);
				seen[tempIndexesLeft[index]] = tempIndexesLeft[index];
			}
			
			//add right if you havent already
			it = seen.find(tempIndexesRight[index]);
			if (it == seen.end()) {  
				mergedResults.push_back(tempIndexesRight[index]);
				seen[tempIndexesRight[index]] = tempIndexesRight[index];
			}
			index++;
		}
		
		
		for (int i = index; i < tempIndexesLeft.size(); i++) {
			if (m->control_pressed) { delete queryRight; delete queryLeft; return refResults; }
			
			//add right if you havent already
			it = seen.find(tempIndexesLeft[i]);
			if (it == seen.end()) {  
				mergedResults.push_back(tempIndexesLeft[i]);
				seen[tempIndexesLeft[i]] = tempIndexesLeft[i];
			}
		}
		
		for (int i = index; i < tempIndexesRight.size(); i++) {
			if (m->control_pressed) { delete queryRight; delete queryLeft; return refResults; }
			
			//add right if you havent already
			it = seen.find(tempIndexesRight[i]);
			if (it == seen.end()) {  
				mergedResults.push_back(tempIndexesRight[i]);
				seen[tempIndexesRight[i]] = tempIndexesRight[i];
			}
		}
		
		for (int i = 0; i < mergedResults.size(); i++) {
			//cout << mergedResults[i]  << '\t' << db[mergedResults[i]]->getName() << endl;	
			if (db[mergedResults[i]]->getName() != q.getName()) { 
				Sequence temp(db[mergedResults[i]]->getName(), db[mergedResults[i]]->getAligned());
				refResults.push_back(temp);
				
			}
		}

		//cout << endl;		
		delete queryRight;
		delete queryLeft;
		
		return refResults;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraSlayer", "getKmerSeqs");
		exit(1);
	}
}
//***************************************************************************************************************


